"""GPU checks (skipped without CUDA): the encoder gives the same planes on the GPU as on the CPU, and the model
trains a step under bf16 autocast."""
import unittest

import torch
import torch.nn as nn

from ..encode import TensorEncoder
from ..env import triton_libcuda
from ..model import PansomaNetV2
from .test_encode import STATS, batch


@unittest.skipUnless(torch.cuda.is_available(), "no CUDA device")
class GpuTest(unittest.TestCase):
    def test_encoder_matches_cpu(self):
        x, blocks = batch(8)
        enc = TensorEncoder(STATS)
        cpu = enc(x, blocks)
        gpu = enc.cuda()(x.cuda(), blocks.cuda()).cpu()
        self.assertTrue(torch.allclose(cpu, gpu, atol=1e-6))

    def test_bf16_training_step(self):
        x, blocks = batch(8)
        model = PansomaNetV2(3, depths=(1, 1, 2, 1), dims=(32, 64, 128, 256), stats=STATS).cuda()
        opt = torch.optim.AdamW(model.parameters(), lr=1e-3)
        y = torch.tensor([0, 1, 2, 0, 1, 2, -1, 0]).cuda()
        loss_fn = nn.CrossEntropyLoss(ignore_index=-1)
        losses = []
        for _ in range(5):
            with torch.autocast("cuda", dtype=torch.bfloat16):
                logits = model(x.cuda(), blocks.cuda())
            loss = loss_fn(logits.float(), y)
            opt.zero_grad()
            loss.backward()
            opt.step()
            losses.append(loss.item())
        self.assertTrue(all(torch.isfinite(torch.tensor(losses))))
        self.assertLess(losses[-1], losses[0])  # it learns the 7 labelled tensors

    def test_channels_last_and_compile_match_eager_in_fp32(self):
        """The speed options change only the kernels: in fp32 without TF32 the compiled channels_last model gives
        the eager model's logits and gradients (compared directly: an AdamW step would turn the sign noise of
        near-zero gradients into full-size updates)."""
        import copy
        triton_libcuda()
        tf32 = (torch.backends.cudnn.allow_tf32, torch.backends.cuda.matmul.allow_tf32)
        torch.backends.cudnn.allow_tf32 = torch.backends.cuda.matmul.allow_tf32 = False
        try:
            x, blocks = batch(8)
            x, blocks = x.cuda(), blocks.cuda()
            y = torch.tensor([0, 1, 2, 0, 1, 2, 0, 1]).cuda()
            torch.manual_seed(0)
            eager = PansomaNetV2(3, depths=(1, 1, 2, 1), dims=(32, 64, 128, 256), stats=STATS, drop_path_rate=0.0).cuda()
            fast = copy.deepcopy(eager).to(memory_format=torch.channels_last)
            compiled = torch.compile(fast)
            with torch.no_grad():
                a, b = eager.eval()(x, blocks), compiled.eval()(x, blocks)
            self.assertLess(float((a - b).abs().max()), 1e-4 * max(1.0, float(a.abs().max())))
            for net, model in ((eager, eager), (compiled, fast)):
                model.train()
                nn.CrossEntropyLoss()(net(x, blocks), y).backward()
            for (name, p), q in zip(eager.named_parameters(), fast.parameters()):
                rel = float((p.grad - q.grad).norm() / (p.grad.norm() + 1e-12))
                self.assertLess(rel, 1e-3, name)
        finally:
            torch.backends.cudnn.allow_tf32, torch.backends.cuda.matmul.allow_tf32 = tf32


if __name__ == "__main__":
    unittest.main()
