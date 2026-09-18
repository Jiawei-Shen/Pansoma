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
