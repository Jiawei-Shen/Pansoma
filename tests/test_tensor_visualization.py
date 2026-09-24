"""Image/metadata regressions; run with a compatible NumPy + Matplotlib environment."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from pangenome_ml_data_generation.tensors import visualization as v


def example():
    x=np.zeros((6,3,4),dtype=np.int16)
    x[:,0]=[[1,6,0,2],[0,-1,0,30],[0,2,2,0],[60,60,0,60],[1,6,0,1],[1,6,6,2]]
    x[:,1]=[[4,3,6,2],[20,25,-1,30],[0,3,2,0],[50,50,50,50],[1,3,6,1],[4,6,6,2]]
    return x,dict(tensor_format_version=v.V2_VERSION,selected_alignments=2,candidate_columns=[1,3],
        coverage=5,alt_count=2,ref_count=1,other_count=2,af=.4)


class VisualizationTest(unittest.TestCase):
    def test_format_requires_v2_identity_and_checks_conflicts(self):
        x,m=example()
        with self.assertRaisesRegex(ValueError,'ambiguous'):
            v.resolve_format(x)
        self.assertEqual(v.resolve_format(x,m),'candidate-v2')
        with self.assertRaises(ValueError):
            v.resolve_format(x,m,'legacy')
        self.assertEqual(v.resolve_format(np.zeros((5,201,100))),'legacy')

    def test_all_channels_trim_and_complete_candidate_range(self):
        x,m=example()
        view,region=v.prepare_candidate_view(x,metadata=m)
        self.assertEqual(view.shape,(6,2,4))
        self.assertEqual(region,[1,3])
        self.assertEqual(v.prepare_candidate_view(x,True,m)[0].shape,(6,3,4))
        with self.assertRaisesRegex(ValueError,'candidate_columns'):
            v.prepare_candidate_view(x,metadata=dict(m,candidate_columns=[1,2]))
        with self.assertRaisesRegex(ValueError,'selected_alignments'):
            v.prepare_candidate_view(x,metadata=dict(m,selected_alignments=1))

    def test_png_preserves_row_zero_qualities_and_per_row_reference(self):
        x,m=example()
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/'v2.png'
            real_close=v.plt.close
            with patch.object(v.plt,'close'):
                self.assertEqual(v.visualize_tensor(x,str(output),'Test',False,None,metadata=m),2)
                fig=v.plt.gcf()
                axes=fig.axes[:6]
                np.testing.assert_array_equal(axes[0].images[0].get_array(),x[0,:2])
                np.testing.assert_array_equal(axes[5].images[0].get_array(),x[5,:2])
                bq=axes[1].images[0].get_array()
                self.assertFalse(bq.mask[0,0])  # BQ=0 on first alignment remains visible.
                self.assertFalse(bq.mask[0,1])  # Gap quality=-1 uses a different color.
                self.assertTrue(bq.mask[0,2])   # Missing coverage is masked white.
                mq=axes[3].images[0].get_array()
                self.assertFalse(mq.mask[0,1])  # MAPQ exists over alignment gaps.
                self.assertEqual(len(axes[0].patches),1)
                self.assertAlmostEqual(axes[0].patches[0].get_width(),2)
                real_close(fig)
            with Image.open(output) as im:
                im.verify()
            self.assertGreater(output.stat().st_size,1000)

    def test_v3_walk_count_panel_and_valid_zero_versus_padding(self):
        x,m=example()
        walks=np.array([[2,2,0,0],[10,10,10,10],[0,0,0,0]],dtype=np.int32)
        seven=np.concatenate([x.astype(np.int32),walks[None]],axis=0)
        m=dict(m,tensor_format_version=v.V3_VERSION)
        self.assertEqual(v.resolve_format(seven,m),"candidate-v3")
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/"v3.png"
            real_close=v.plt.close
            with patch.object(v.plt,"close"):
                v.visualize_tensor(seven,str(output),"V3",False,None,metadata=m)
                fig=v.plt.gcf()
                values=fig.axes[6].images[0].get_array()
                self.assertEqual(int(values[0,0]),2)
                self.assertTrue(values.mask[0,2])
                self.assertFalse(values.mask[0,3])  # Known node with zero W count.
                self.assertEqual(int(values[0,3]),0)
                real_close(fig)
            with Image.open(output) as image:
                image.verify()

    def test_v4_explicit_format_and_gbwt_legend(self):
        x,m=example()
        x=np.concatenate([x.astype(np.int32),np.full_like(x[:1],86)],axis=0)
        x[6][x[4]==0]=0
        m=dict(m,tensor_format_version=v.V4_VERSION)
        self.assertEqual(v.resolve_format(x,m),"candidate-v4")
        with self.assertRaisesRegex(ValueError,"conflicts"):
            v.resolve_format(x,m,"candidate-v3")
        with tempfile.TemporaryDirectory() as d:
            real_close=v.plt.close
            with patch.object(v.plt,"close"):
                v.visualize_tensor(x,str(Path(d)/"v4.png"),"V4",False,None,tensor_format="candidate-v4")
                fig=v.plt.gcf()
                labels=[t.get_text() for ax in fig.axes for t in ax.texts]
                self.assertIn("candidate-v4",fig._suptitle.get_text())
                self.assertIn("Distinct GBWT\npaths",labels)
                self.assertNotIn("Distinct GFA\nW records",labels)
                real_close(fig)

    def test_int8_v4_renders_scaled_counts(self):
        x,m=example()
        x=np.concatenate([x.astype(np.int32),np.full_like(x[:1],86)],axis=0)
        x = np.clip(x, -1, 127).astype(np.int8)
        x[6] = 82
        x[6][x[4]==0]=0
        m=dict(m,tensor_format_version=v.V4_VERSION, tensor_storage_version="int8-count-div4-v1")
        self.assertEqual(v.resolve_format(x,m),"candidate-v4")
        with self.assertRaisesRegex(ValueError,"conflicts"):
            v.resolve_format(x,m,"candidate-v3")
        with tempfile.TemporaryDirectory() as d:
            real_close=v.plt.close
            with patch.object(v.plt,"close"):
                v.visualize_tensor(x,str(Path(d)/"v4.png"),"V4",False,None,tensor_format="candidate-v4", metadata=m)
                fig=v.plt.gcf()
                self.assertEqual(int(fig.axes[6].images[0].get_array().max()), 82)
                labels=[t.get_text() for ax in fig.axes for t in ax.texts]
                self.assertIn("candidate-v4",fig._suptitle.get_text())
                self.assertIn("Distinct GBWT\npaths\n// 4 (cap 127)",labels)
                self.assertNotIn("Distinct GFA\nW records",labels)
                real_close(fig)

    def test_v5_alt_stripe_log_counts_and_strand_panels(self):
        x,m=example()
        x=x.astype(np.int8)
        evidence=x[0]!=0
        stripe=np.zeros_like(x[0]); stripe[:, 1:3]=[[4, 1]]*3; stripe[~evidence]=0   # candidate INS ">TA" at columns 1-2
        counts=np.where(evidence, 91, 0).astype(np.int8)                             # 90 paths -> 91 on the log scale
        strand=np.where(evidence, [[1],[2],[0]], 0).astype(np.int8)
        v5=np.concatenate([x[:2], stripe[None], x[3:6], counts[None], strand[None]], axis=0)
        v5[2][~evidence]=0
        m=dict(m, tensor_format_version=v.V5_VERSION, tensor_storage_version=v.V5_LOG_STORAGE, candidate_columns=[1,3])
        self.assertEqual(v.resolve_format(v5,m),"candidate-v5")
        self.assertEqual(v.resolve_format(np.zeros((8,3,4),dtype=np.int8),tensor_format="candidate-v5"),"candidate-v5")
        with self.assertRaisesRegex(ValueError,"8 channels"):
            v.resolve_format(x,dict(m))
        view,region=v.prepare_candidate_view(v5,metadata=m)
        self.assertEqual((view.shape,region),((8,2,4),[1,3]))
        # A stripe that ends early (no evidence in a trailing candidate column) is still accepted.
        short=v5.copy(); short[2,:,2]=0
        self.assertEqual(v.prepare_candidate_view(short,metadata=m)[1],[1,3])
        with self.assertRaisesRegex(ValueError,"candidate_columns"):
            v.prepare_candidate_view(v5,metadata=dict(m,candidate_columns=[2,3]))
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/"v5.png"
            real_close=v.plt.close
            with patch.object(v.plt,"close"):
                self.assertEqual(v.visualize_tensor(v5,str(output),"V5",False,None,metadata=m),2)
                fig=v.plt.gcf()
                panels=[ax for ax in fig.axes if ax.images]
                self.assertEqual(len(panels),8)
                np.testing.assert_array_equal(panels[2].images[0].get_array(),v5[2,:2])   # ALT stripe in base palette
                np.testing.assert_array_equal(panels[7].images[0].get_array(),v5[7,:2])   # strand
                self.assertEqual(int(panels[6].images[0].get_array().max()),91)
                labels=[t.get_text() for ax in fig.axes for t in ax.texts]
                self.assertIn("Distinct GBWT\npaths\nlog scale; ticks\nshow counts",labels)
                titles=[ax.get_title(loc="left") for ax in panels]
                self.assertIn("3  Candidate ALT (same for every row)",titles)
                self.assertIn("8  Read strand (vs. candidate node forward)",titles)
                self.assertIn("candidate-v5",fig._suptitle.get_text())
                real_close(fig)
            with Image.open(output) as image:
                image.verify()

    def test_v6_site_allele_blocks_and_caption(self):
        x,m=example()
        x=x.astype(np.int8)
        evidence=x[0]!=0
        stripe=np.zeros_like(x[0]); stripe[0,1:3]=[4,1]; stripe[1,1:3]=[6,6]; stripe[~evidence]=0  # row 0 A1 ">TA", row 1 REF
        counts=np.where(evidence,91,0).astype(np.int8)
        strand=np.where(evidence,[[1],[2],[0]],0).astype(np.int8)
        v6=np.concatenate([x[:2],stripe[None],x[3:6],counts[None],strand[None]],axis=0)
        allele=dict(label="A1",event_type="INS",ref="",alt="TA",af=.4,alt_count=2,coverage=5)
        m=dict(m,tensor_format_version=v.V6_VERSION,tensor_storage_version=v.V6_LINEAR_STORAGE,candidate_columns=[1,3],
               alleles=[allele],site_coverage=6,site_counts=dict(A1=2,REF=1,OTHER=3),
               row_groups=[dict(start_row=0,end_row=1,allele="A1"),dict(start_row=1,end_row=2,allele="REF")])
        self.assertEqual(v.resolve_format(v6,m),"candidate-v6")
        self.assertEqual(v.prepare_candidate_view(v6,metadata=m)[1],[1,3])
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/"v6.png"
            real_close=v.plt.close
            with patch.object(v.plt,"close"):
                self.assertEqual(v.visualize_tensor(v6,str(output),"V6",False,None,metadata=m),2)
                fig=v.plt.gcf()
                panels=[ax for ax in fig.axes if ax.images]
                np.testing.assert_array_equal(panels[2].images[0].get_array(),v6[2,:2])   # per-row site allele
                titles=[ax.get_title(loc="left") for ax in panels]
                self.assertIn("3  Site allele carried by each read (A1, A2, ... or REF; blank = OTHER)",titles)
                texts=[t.get_text() for ax in panels for t in ax.texts]
                self.assertEqual(texts.count("A1"),2)  # block labels on panels 1 and 3
                self.assertEqual(texts.count("REF"),2)
                labels=[t.get_text() for ax in fig.axes for t in ax.texts]
                self.assertIn("Distinct GBWT\npaths\nexact <= 100,\nlog2 above",labels)
                caption=fig._suptitle.get_text()
                self.assertIn("A1 INS ->TA AF 0.400 (2/5)",caption)
                self.assertIn("Site coverage 6 (A1 2, REF 1, OTHER 3)",caption)
                real_close(fig)
            with Image.open(output) as image:
                image.verify()

    def test_legacy_still_renders_five_channels(self):
        x=np.zeros((5,3,4),dtype=np.int8)
        x[0,:2]=20
        x[1,1]=30
        x[2,1,2]=5
        x[3,1]=60
        x[4,1]=10
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/'legacy.png'
            self.assertEqual(v.visualize_tensor(x,str(output),'Legacy',False,2),2)
            with Image.open(output) as im:
                im.verify()

    def test_metadata_auto_loading_and_custom_shard_index(self):
        _,m=example()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'manifest.json').write_text(json.dumps({'tensor_format_version':v.V2_VERSION}))
            record=dict(m,shard_index=3,index_within_shard=0)
            (p/'variant_summary.ndjson').write_text(json.dumps(record)+'\n')
            manifest,rows=v.load_metadata(str(p/'shard_00003_data.npy'))
            self.assertEqual(manifest['tensor_format_version'],v.V2_VERSION)
            self.assertEqual(rows[0]['candidate_columns'],[1,3])
            self.assertEqual(v.load_metadata(str(p/'custom.npy'),shard_index=3)[1],rows)
            (p/'variant_summary.ndjson').write_text((json.dumps(record)+'\n')*2)
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                v.load_metadata(str(p/'custom.npy'),shard_index=3)


if __name__=='__main__':
    unittest.main()
