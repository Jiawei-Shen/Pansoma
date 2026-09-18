#!/usr/bin/env python3
"""Compare saved legacy/candidate (v2 or v3) summaries; independently audit SNP counts in source GAM."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from indexed_gam_pipeline.gam_reader import IndexedGam
from indexed_gam_pipeline.candidates import rc
from indexed_gam_pipeline.run import node_records


def legacy_identity(meta):
    kind = {"X": "SNP", "I": "INS", "D": "DEL"}[meta["v_type"]]
    # Legacy insertion position is the preceding base; v_alt is unanchored,
    # although variant_key includes an anchored REF/ALT spelling.
    pos = meta["v_pos"] + (kind == "INS")
    ref = "" if kind == "INS" else meta["v_ref"]
    alt = "" if kind == "DEL" else meta["v_alt"]
    return (int(meta["node_id"]), pos, kind, ref, alt)


def v2_identity(meta):
    return (int(meta["node_id"]), meta["start"], meta["event_type"], meta["ref"], meta["alt"])


def audit_snp(alignment, candidate, node_length, min_bq):
    """Directly walk raw edit/read cursors, independently of candidates.overlap."""
    read_cursor = 0
    hits = []
    for mi, mapping in enumerate(alignment.path.mapping):
        cursor = mapping.position.offset
        reverse = mapping.position.is_reverse
        target = node_length-1-candidate["start"] if reverse else candidate["start"]
        for ei, edit in enumerate(mapping.edit):
            f, t = edit.from_length, edit.to_length
            if mapping.position.node_id == candidate["node_id"] and cursor <= target < cursor+f:
                base, bq, op = "-", None, "D" if t == 0 else "complex"
                if f == t:
                    index = read_cursor+target-cursor
                    base = alignment.sequence[index].upper()
                    if reverse:
                        base = rc(base)
                    bq = alignment.quality[index] if alignment.quality else -1
                    op = "X" if edit.sequence else "M"
                old = "alt" if base == candidate["alt"] else "ref" if base == candidate["ref"] else "other"
                new = ("alt" if op == "X" and base == candidate["alt"] and bq >= min_bq
                       else "ref" if op in ("M", "X") and base == candidate["ref"] else "other")
                hits.append(dict(legacy=old, candidate_v2=new, forward_base=base, base_quality=bq,
                    operation=op, mapping_index=mi, edit_index=ei, mapping_reverse=reverse,
                    edit=dict(from_length=f,to_length=t,sequence=edit.sequence)))
            cursor += f
            read_cursor += t
    if not hits:
        return None
    hit = min(hits, key=lambda h: ({"alt":0,"ref":1,"other":2}[h["candidate_v2"]],h["mapping_index"]))
    return dict(hit, read_name=alignment.name, mapping_quality=alignment.mapping_quality,
                record_sha256=hashlib.sha256(alignment.SerializeToString(deterministic=True)).hexdigest(),
                visits=len(hits))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--legacy',required=True,help='Legacy summary NDJSON')
    parser.add_argument('--candidate', '--v2', dest='v2', required=True, help='Candidate v2/v3 summary NDJSON (--v2 retained as an alias)')
    parser.add_argument('--gam',required=True)
    parser.add_argument('--index')
    parser.add_argument('--node-sqlite',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    old={legacy_identity(m):m for m in map(json.loads,Path(args.legacy).read_text().splitlines())}
    new={v2_identity(m):m for m in map(json.loads,Path(args.v2).read_text().splitlines())}
    nodes={k[0] for k in new}
    records=node_records(argparse.Namespace(node_json=None,gfa=None,node_sqlite=args.node_sqlite),nodes)
    snps={k:[] for k in new if k[2]=='SNP'}
    for a in IndexedGam(args.gam,args.index).fetch(nodes):
        for key,hits in snps.items():
            meta=new[key]
            if a.mapping_quality <= meta['parameters']['min_mapq']:
                continue
            hit=audit_snp(a,meta,len(records[key[0]]['sequence']),meta['parameters']['min_allele_bq'])
            if hit:
                hits.append(hit)
    comparisons=[]
    for key in sorted(old.keys() | new.keys()):
        a,b=old.get(key),new.get(key)
        item=dict(identity=key,matched=a is not None and b is not None)
        if a:
            item['legacy']=dict(coverage=a['coverage_at_locus'],alt=a['alt_allele_count'],
                ref=a['ref_allele_count_at_locus'],other=a['other_allele_count_at_locus'],
                af=a['alt_allele_frequency'],variant_key=a['variant_key'])
        if b:
            item['candidate_v2']=dict(coverage=b['coverage'],alt=b['alt_count'],ref=b['ref_count'],
                other=b['other_count'],af=b['af'],candidate_id=b['candidate_id'])
        if a and b:
            item['count_deltas_v2_minus_legacy']={field:item['candidate_v2'][field]-item['legacy'][field]
                                               for field in ('coverage','alt','ref','other')}
            item['af_equal_at_legacy_precision']=round(b['af'],4)==a['alt_allele_frequency']
        if key in snps:
            hits=snps[key]
            before=Counter(h['legacy'] for h in hits)
            after=Counter(h['candidate_v2'] for h in hits)
            audit=dict(coverage=len(hits),legacy=dict(before),candidate_v2=dict(after),
                       changed_observations=[h for h in hits if h['legacy']!=h['candidate_v2']])
            audit['v2_matches_raw_edits']=(len(hits)==b['coverage'] and all(after[c]==b[c+'_count'] for c in ('alt','ref','other')))
            audit['legacy_matches_raw_bases']=(a is not None and len(hits)==a['coverage_at_locus'] and
                all(before[c]==item['legacy'][c] for c in ('alt','ref','other')))
            item['raw_snp_audit']=audit
        comparisons.append(item)
    report=dict(legacy_summary=args.legacy,v2_summary=args.v2,comparisons=comparisons,
        raw_snp_audits_passed=all(c['raw_snp_audit']['v2_matches_raw_edits'] and c['raw_snp_audit']['legacy_matches_raw_bases']
                                for c in comparisons if 'raw_snp_audit' in c),
        note='Legacy minimum BQ filters mean ALT BQ; v2 applies minimum BQ per ALT observation. Legacy AF is rounded to four decimals.')
    versions = {m['tensor_format_version'] for m in new.values()}
    if len(versions) != 1:
        raise ValueError('Comparison requires one candidate tensor format')
    report['tensor_format_version'] = next(iter(versions))
    if report['tensor_format_version'] in ('indexed-gam-candidate-v3', 'indexed-gam-candidate-v4'):
        target_version = report['tensor_format_version'].rsplit('-', 1)[-1]
        # Keep v2 reports backward compatible while labeling the current format accurately.
        def relabel(value):
            if isinstance(value, dict):
                return {k.replace('v2', target_version): relabel(v) for k, v in value.items()}
            if isinstance(value, list):
                return [relabel(v) for v in value]
            return value
        report = relabel(report)
        report['note'] = report['note'].replace('v2', target_version)
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if not report['raw_snp_audits_passed']:
        raise SystemExit('Raw SNP audit failed')


if __name__=='__main__':
    main()
