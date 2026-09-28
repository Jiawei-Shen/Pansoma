"""PansomaNetV2 = TensorEncoder (encode.py) -> 1x1 front -> ConvNeXt-CBAM backbone -> 3 classes.

The backbone (LayerNorm, GRN, CBAM, ConvNeXtBlock, DownsampleLayer, ConvNeXtCBAM) is the one of
machine_learning/pansoma_net/mynet.py (ConvNeXtCBAMClassifier), copied unchanged except that it no longer
prints its configuration, so v2 does not depend on v1's files. The 1x1 front (two 1x1 convolutions with
GELU by default) combines the 36 planes of each cell before the 4 x 4, stride 4 stem mixes cells.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from timm.layers import DropPath, trunc_normal_

from .encode import N_PLANES, TensorEncoder


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
    def __init__(self, dim):
        super().__init__()
        self.gamma = nn.Parameter(torch.zeros(1))
        self.beta = nn.Parameter(torch.zeros(1))

    def forward(self, x):
        gx = torch.norm(x, dim=(2, 3), keepdim=True)
        nx = gx / (gx.mean(dim=1, keepdim=True) + 1e-6)
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
    def __init__(self, in_channels, out_channels, drop_path):
        super().__init__()
        self.dwconv = nn.Conv2d(in_channels, in_channels, kernel_size=7, padding=3, groups=in_channels)
        self.norm = LayerNorm(in_channels, eps=1e-6)
        self.pwconv1 = nn.Linear(in_channels, 4 * in_channels)
        self.act = nn.GELU()
        self.grn = GRN(4 * in_channels)
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
        return self.act(self.proj(shortcut) + x)


class DownsampleLayer(nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.downsample = nn.Sequential(LayerNorm(in_channels, eps=1e-6, data_format="channels_first"),
                                        nn.Conv2d(in_channels, out_channels, kernel_size=2, stride=2), nn.GELU())

    def forward(self, x):
        return self.downsample(x)


class ConvNeXtCBAM(nn.Module):
    """mynet.ConvNeXtCBAMClassifier: 4 x 4 stride 4 stem, four stages, three 2 x 2 downsamplings, pooled head."""

    def __init__(self, in_channels, class_num, depths=(3, 3, 27, 3), dims=(192, 384, 768, 1536), drop_path_rate=0.1):
        super().__init__()
        self.stem = nn.Sequential(nn.Conv2d(in_channels, dims[0], kernel_size=4, stride=4),
                                  LayerNorm(dims[0], eps=1e-6, data_format="channels_first"), nn.GELU())
        rates = [x.item() for x in torch.linspace(0, drop_path_rate, sum(depths))]
        self.stages, self.downsample_layers, k = nn.ModuleList(), nn.ModuleList(), 0
        for i in range(4):
            self.stages.append(nn.Sequential(*[ConvNeXtBlock(dims[i], dims[i], rates[k + j]) for j in range(depths[i])]))
            k += depths[i]
        for i in range(3):
            self.downsample_layers.append(DownsampleLayer(dims[i], dims[i + 1]))
        self.pool = nn.AdaptiveAvgPool2d((1, 1))
        self.head = nn.Linear(dims[-1], class_num)
        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, (nn.Conv2d, nn.Linear)):
            trunc_normal_(m.weight, std=.02)
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)

    def forward(self, x):
        x = self.stem(x)
        for i in range(4):
            x = self.stages[i](x)
            if i < 3:
                x = self.downsample_layers[i](x)
        return self.head(torch.flatten(self.pool(x), 1))


class PansomaNetV2(nn.Module):
    def __init__(self, num_classes=3, depths=(3, 3, 27, 3), dims=(192, 384, 768, 1536), front=(64, 64),
                 drop_path_rate=0.1, stats=None):
        super().__init__()
        self.config = dict(num_classes=num_classes, depths=list(depths), dims=list(dims), front=list(front),
                           drop_path_rate=drop_path_rate)
        self.encoder = TensorEncoder(stats)
        layers, c = [], N_PLANES
        for width in front:
            layers += [nn.Conv2d(c, width, kernel_size=1), nn.GELU()]
            c = width
        self.front = nn.Sequential(*layers)
        self.backbone = ConvNeXtCBAM(c, num_classes, depths, dims, drop_path_rate)
        for m in self.front.modules():
            if isinstance(m, nn.Conv2d):
                trunc_normal_(m.weight, std=.02)
                nn.init.constant_(m.bias, 0)

    def forward(self, x, blocks):
        """x: (B, 8, 200, 101) int8 tensors; blocks: (B, 4) row-block ends -> logits (B, num_classes)."""
        return self.backbone(self.front(self.encoder(x, blocks)))

    @classmethod
    def from_checkpoint(cls, checkpoint, map_location="cpu"):
        if not isinstance(checkpoint, dict):
            checkpoint = torch.load(checkpoint, map_location=map_location, weights_only=False)
        if checkpoint.get("format") != "pansoma_net_v2":
            raise ValueError("not a pansoma_net_v2 checkpoint")
        model = cls(**checkpoint["config"])
        model.load_state_dict(checkpoint["model_state_dict"])
        return model
