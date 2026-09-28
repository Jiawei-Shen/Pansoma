"""GPU checks (skipped without CUDA): the encoder gives the same planes on the GPU as on the CPU, and the model
trains a step under bf16 autocast."""
import unittest

import torch
import torch.nn as nn

from ..encode import TensorEncoder
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


if __name__ == "__main__":
    unittest.main()
