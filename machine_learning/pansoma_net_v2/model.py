"""PansomaNetV2 = TensorEncoder (encode.py) -> 1x1 front -> ConvNeXt-CBAM backbone -> 3 classes.

The backbone (LayerNorm, GRN, CBAM, ConvNeXtBlock, DownsampleLayer, ConvNeXtCBAM) comes from
machine_learning/pansoma_net/mynet.py (ConvNeXtCBAMClassifier), so v2 does not depend on v1's files. The 1x1
front (two 1x1 convolutions with GELU by default) combines the 36 planes of each cell before the 4 x 4, stride 4
stem mixes cells.

config "block" chooses the backbone's blocks:
- "v2" (the default): ConvNeXt(-V2) as published. GRN takes each channel's L2 norm over the spatial positions
  and divides by its mean over the channels, with per-channel gamma / beta; a block returns shortcut + branch
  (an identity path through every stage); the stem and the downsampling layers have no GELU.
- "v1": mynet's blocks unchanged, kept to load and compare earlier checkpoints (a config without "block").
  Its GRN normalizes each read row by its norm over (positions, channels) with one scalar gamma / beta, and a
  GELU follows every residual sum, the stem and every downsampling.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from timm.layers import DropPath, trunc_normal_

from .encode import TensorEncoder


class LayerNorm(nn.Module):
    """LayerNorm for channels_last (B, H, W, C) or channels_first (B, C, H, W) inputs."""

    def __init__(self, normalized_shape, eps=1e-6, data_format="channels_last"):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(normalized_shape))
        self.bias = nn.Parameter(torch.zeros(normalized_shape))
        self.eps = eps
        self.data_format = data_format
        if self.data_format not in ["channels_last", "channels_first"]:
            raise NotImplementedError
        self.normalized_shape = (normalized_shape,)

    def forward(self, x):
        if self.data_format == "channels_last":
            return F.layer_norm(x, self.normalized_shape, self.weight, self.bias, self.eps)
        u = x.mean(1, keepdim=True)
        s = (x - u).pow(2).mean(1, keepdim=True)
        x = (x - u) / torch.sqrt(s + self.eps)
        return self.weight[:, None, None] * x + self.bias[:, None, None]


class GRN(nn.Module):
    """mynet's GRN (block "v1"): on a channels_last (B, H, W, C) input it takes the norm over (W, C), i.e. per
    read row, with scalar gamma / beta."""

    def __init__(self, dim):
        super().__init__()
        self.gamma = nn.Parameter(torch.zeros(1))
        self.beta = nn.Parameter(torch.zeros(1))

    def forward(self, x):
        gx = torch.norm(x, dim=(2, 3), keepdim=True)
        nx = gx / (gx.mean(dim=1, keepdim=True) + 1e-6)
        return self.gamma * (x * nx) + self.beta + x


class GRNv2(nn.Module):
    """Global Response Normalization of ConvNeXt V2 (channels_last (B, H, W, C)): per channel the L2 norm over
    the positions, divided by its mean over the channels; per-channel gamma / beta, zero at the start (identity)."""

    def __init__(self, dim):
        super().__init__()
        self.gamma = nn.Parameter(torch.zeros(1, 1, 1, dim))
        self.beta = nn.Parameter(torch.zeros(1, 1, 1, dim))

    def forward(self, x):
        gx = torch.norm(x, p=2, dim=(1, 2), keepdim=True)
        nx = gx / (gx.mean(dim=-1, keepdim=True) + 1e-6)
        return self.gamma * (x * nx) + self.beta + x


class ChannelAttention(nn.Module):
    def __init__(self, in_planes, ratio=16):
        super().__init__()
        self.avg_pool = nn.AdaptiveAvgPool2d(1)
        self.max_pool = nn.AdaptiveMaxPool2d(1)
        self.fc = nn.Sequential(nn.Conv2d(in_planes, in_planes // ratio, 1), nn.ReLU(),
                                nn.Conv2d(in_planes // ratio, in_planes, 1))
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        return x * self.sigmoid(self.fc(self.avg_pool(x)) + self.fc(self.max_pool(x)))


class SpatialAttention(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = nn.Conv2d(2, 1, kernel_size=7, padding=3)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        attn = torch.cat([torch.mean(x, dim=1, keepdim=True), torch.max(x, dim=1, keepdim=True)[0]], dim=1)
        return x * self.sigmoid(self.conv(attn))


class CBAM(nn.Module):
    def __init__(self, channels):
        super().__init__()
        self.ca = ChannelAttention(channels)
        self.sa = SpatialAttention()

    def forward(self, x):
        return self.sa(self.ca(x))


class ConvNeXtBlock(nn.Module):
    def __init__(self, in_channels, out_channels, drop_path, block="v1"):
        super().__init__()
        self.v2 = block == "v2"
        self.dwconv = nn.Conv2d(in_channels, in_channels, kernel_size=7, padding=3, groups=in_channels)
        self.norm = LayerNorm(in_channels, eps=1e-6)
        self.pwconv1 = nn.Linear(in_channels, 4 * in_channels)
        self.act = nn.GELU()
        self.grn = (GRNv2 if self.v2 else GRN)(4 * in_channels)
        self.pwconv2 = nn.Linear(4 * in_channels, out_channels)
        self.cbam = CBAM(out_channels)
        self.drop_path = DropPath(drop_path) if drop_path > 0. else nn.Identity()
        self.proj = nn.Identity() if in_channels == out_channels else nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=1), nn.GELU())

    def forward(self, x):
        shortcut = x
        x = self.dwconv(x).permute(0, 2, 3, 1)
        x = self.pwconv2(self.grn(self.act(self.pwconv1(self.norm(x))))).permute(0, 3, 1, 2)
        x = self.drop_path(self.cbam(x))
        if self.v2:
            return self.proj(shortcut) + x
        return self.act(self.proj(shortcut) + x)


class DownsampleLayer(nn.Module):
    def __init__(self, in_channels, out_channels, block="v1"):
        super().__init__()
        layers = [LayerNorm(in_channels, eps=1e-6, data_format="channels_first"),
                  nn.Conv2d(in_channels, out_channels, kernel_size=2, stride=2)]
        self.downsample = nn.Sequential(*layers, *([] if block == "v2" else [nn.GELU()]))

    def forward(self, x):
        return self.downsample(x)


class ConvNeXtCBAM(nn.Module):
    """mynet.ConvNeXtCBAMClassifier: 4 x 4 stride 4 stem, four stages, three 2 x 2 downsamplings, pooled head."""

    def __init__(self, in_channels, class_num, depths=(3, 3, 27, 3), dims=(192, 384, 768, 1536), drop_path_rate=0.1,
                 block="v1"):
        super().__init__()
        if block not in ("v1", "v2"):
            raise ValueError(f"block {block!r}: v1 or v2")
        self.stem = nn.Sequential(nn.Conv2d(in_channels, dims[0], kernel_size=4, stride=4),
                                  LayerNorm(dims[0], eps=1e-6, data_format="channels_first"),
                                  *([] if block == "v2" else [nn.GELU()]))
        rates = [x.item() for x in torch.linspace(0, drop_path_rate, sum(depths))]
        self.stages, self.downsample_layers, k = nn.ModuleList(), nn.ModuleList(), 0
        for i in range(4):
            self.stages.append(nn.Sequential(*[ConvNeXtBlock(dims[i], dims[i], rates[k + j], block)
                                               for j in range(depths[i])]))
            k += depths[i]
        for i in range(3):
            self.downsample_layers.append(DownsampleLayer(dims[i], dims[i + 1], block))
        self.pool = nn.AdaptiveAvgPool2d((1, 1))
        self.head = nn.Linear(dims[-1], class_num)
        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, (nn.Conv2d, nn.Linear)):
            trunc_normal_(m.weight, std=.02)
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)

    def forward_features(self, x):
        x = self.stem(x)
        for i in range(4):
            x = self.stages[i](x)
            if i < 3:
                x = self.downsample_layers[i](x)
        return torch.flatten(self.pool(x), 1)

    def forward(self, x):
        return self.head(self.forward_features(x))


def full_config(config):
    """A checkpoint's config with the defaults of keys added later (checkpoints before "block" are v1)."""
    return {"scalars": 0, "block": "v1", "drop_planes": [], **config}


def no_decay(name, param):
    """Parameters AdamW should not decay: biases, norm weights (1-D) and the GRN gamma / beta."""
    return param.ndim <= 1 or ".grn." in name


class PansomaNetV2(nn.Module):
    """scalars > 0: the site's scalars (data.SCALARS), z-scored with statistics fitted on the training tensors
    (buffers, like the encoder's), go through a Linear(scalars, 64) + GELU and join the pooled features before the
    head."""

    def __init__(self, num_classes=3, depths=(3, 3, 27, 3), dims=(192, 384, 768, 1536), front=(64, 64),
                 drop_path_rate=0.1, stats=None, scalars=0, block="v2", drop_planes=()):
        super().__init__()
        self.config = dict(num_classes=num_classes, depths=list(depths), dims=list(dims), front=list(front),
                           drop_path_rate=drop_path_rate, scalars=scalars, block=block, drop_planes=list(drop_planes))
        self.encoder = TensorEncoder(stats, drop_planes)
        layers, c = [], self.encoder.n_planes
        for width in front:
            layers += [nn.Conv2d(c, width, kernel_size=1), nn.GELU()]
            c = width
        self.front = nn.Sequential(*layers)
        self.backbone = ConvNeXtCBAM(c, num_classes, depths, dims, drop_path_rate, block)
        if scalars:
            self.register_buffer("scalar_mean", torch.zeros(scalars))
            self.register_buffer("scalar_std", torch.ones(scalars))
            self.scalar_mlp = nn.Sequential(nn.Linear(scalars, 64), nn.GELU())
            self.backbone.head = nn.Linear(dims[-1] + 64, num_classes)
        for m in list(self.front.modules()) + list(self.scalar_mlp.modules() if scalars else []) + [self.backbone.head]:
            if isinstance(m, (nn.Conv2d, nn.Linear)):
                trunc_normal_(m.weight, std=.02)
                nn.init.constant_(m.bias, 0)

    def set_scalar_stats(self, mean, std):
        self.scalar_mean.copy_(torch.as_tensor(mean, dtype=torch.float32))
        self.scalar_std.copy_(torch.clamp(torch.as_tensor(std, dtype=torch.float32), min=1e-6))

    def forward(self, x, blocks, scalars=None):
        """x: (B, 8, 200, 101) int8 tensors; blocks: (B, 4) row-block ends; scalars: (B, S) -> logits (B, classes)."""
        features = self.backbone.forward_features(self.front(self.encoder(x, blocks)))
        if self.config["scalars"]:
            z = (scalars.to(features.device, torch.float32) - self.scalar_mean) / self.scalar_std
            features = torch.cat([features, self.scalar_mlp(z).to(features.dtype)], 1)
        return self.backbone.head(features)

    @classmethod
    def from_checkpoint(cls, checkpoint, map_location="cpu"):
        if not isinstance(checkpoint, dict):
            checkpoint = torch.load(checkpoint, map_location=map_location, weights_only=False)
        if checkpoint.get("format") != "pansoma_net_v2":
            raise ValueError("not a pansoma_net_v2 checkpoint")
        model = cls(**full_config(checkpoint["config"]))
        model.load_state_dict(checkpoint["model_state_dict"])
        return model
