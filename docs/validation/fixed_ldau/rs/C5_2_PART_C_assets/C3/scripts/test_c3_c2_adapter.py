import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import c3_c2_adapter as c3


def write(path: Path, data: bytes) -> str:
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


class C3AdapterSmoke(unittest.TestCase):
    def test_preregister_accepts_baseline_writer_shape_and_rejects_wrong_bindings(self):
        # Match the nested `capture` object written by baseline mode, then make
        # each identity mismatch fail closed before any replay cases are run.
        capture={'csv_sha256':'capture-csv-sha','manifest_sha256':'capture-manifest-sha'}
        report={'schema':'rslmto-c5-c3-live-baseline-v1',
                'status':'PASS_LIVE_STATE_BASELINE','capture':capture,
                'driver_sha256':'driver-sha',
                'baseline_replay_integrity_checks':{'pass':True}}
        c3.validate_live_state_baseline_binding(report,capture,'driver-sha')

        wrong_csv={**report,'capture':{**capture,'csv_sha256':'different-csv'}}
        with self.assertRaisesRegex(ValueError,'not bound to this capture/driver'):
            c3.validate_live_state_baseline_binding(wrong_csv,capture,'driver-sha')
        wrong_manifest={**report,'capture':{**capture,'manifest_sha256':'different-manifest'}}
        with self.assertRaisesRegex(ValueError,'not bound to this capture/driver'):
            c3.validate_live_state_baseline_binding(wrong_manifest,capture,'driver-sha')
        with self.assertRaisesRegex(ValueError,'not bound to this capture/driver'):
            c3.validate_live_state_baseline_binding(report,capture,'different-driver')
        with self.assertRaisesRegex(ValueError,'missing nested capture binding'):
            c3.validate_live_state_baseline_binding(
                {'status':'PASS_LIVE_STATE_BASELINE','capture_csv_sha256':'capture-csv-sha',
                 'manifest_sha256':'capture-manifest-sha','driver_sha256':'driver-sha'},
                capture,'driver-sha')

    def test_requires_accepted_c2_chebyshev_endpoint(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            deck=root/'input.nml'; seed=root/'Fe_out.nml'; binary=root/'rslmto.x'; source=root/'self.f90'
            deck_text=(b"&self\nuse_kspace=.false.\n/\n&control\ncalctype='B'\nnsp=1\n"
                       b"recur='chebyshev'\n/\n&hamiltonian\nhubbard_radial_moments='local'\n/")
            dsha=write(deck,deck_text); ssha=write(seed,b'seed'); bsha=write(binary,b'binary'); fsha=write(source,b'source')
            conv={'outer_diff':2e-8,'last20_d_occupation_by_spin':{'1':[1.0]*20,'2':[2.0]*20},
                  'last20_sphere_moment_muB':[3.0]*20}
            lock={'schema':c3.LOCK_SCHEMA,'status':c3.LOCK_STATUS,
                  'system':{'backend':'RS','calctype':'B','recurrence':'chebyshev','radial_moments':'local',
                            'nsp':1,'lmax':2,'u_ev':2.0,'j_ev':0.0},
                  'acceptance':{'c2_gate_status':'NOT_EVALUATED','c1_rs_status':'ACCEPTED','convergence':conv},
                  'accepted_endpoint':{'run_id':'synthetic','accepted_iteration':21,
                       'deck_path':str(deck),'deck_sha256':dsha,'checkpoint_path':str(seed),'checkpoint_sha256':ssha},
                  'executable':{'binary_path':str(binary),'binary_sha256':bsha,
                                'source_path':str(source),'source_sha256':fsha}}
            lock_path=root/'lock.json'; lock_path.write_text(json.dumps(lock))
            binding=c3.validate_c2_lock(lock,lock_path)
            self.assertEqual(binding['system']['recurrence'],'chebyshev')
            self.assertEqual(binding['acceptance']['block_chebyshev_comparison_status'],'NOT_EVALUATED')
            self.assertNotIn('co_u_positive_admitted',binding)
            lock['system']['recurrence']='block'
            with self.assertRaisesRegex(ValueError,'Chebyshev'):
                c3.validate_c2_lock(lock,lock_path)
            lock['system']['recurrence']='chebyshev'
            invalid=deck_text.replace(b"&self\nuse_kspace=.false.\n/",b"&self\n/\n&control\nuse_kspace=.false.")
            dsha=write(deck,invalid); lock['accepted_endpoint']['deck_sha256']=dsha
            with self.assertRaisesRegex(ValueError,'unknown members=.*use_kspace'):
                c3.validate_c2_lock(lock,lock_path)

    def test_full_green_capture_requires_every_matrix_pair(self):
        rows=[]
        labels=[(0,0,0),(0,0,1),(1,-1,0),(1,0,0),(1,1,0),
                (1,-1,1),(1,0,1),(1,1,1)]
        for i,(ell,m,spin) in enumerate(labels):
            rows.append({'kind':'basis_map','stage':'RS_GREEN','site':'0','row':str(i),
                         'l':str(ell),'col':str(m),'spin':str(spin)})
        for energy in range(3):
            for i in range(8):
                for j in range(8):
                    rows.append({'kind':'green_function','stage':'RS_GREEN','site':'0','k':str(energy),
                                 'row':str(i),'col':str(j)})
        shape=c3.validate_full_green_rows(rows)
        self.assertEqual(shape['energy_grid_points_by_site'],{'0':3})
        del rows[-1]
        with self.assertRaisesRegex(ValueError,'incomplete Green matrix'):
            c3.validate_full_green_rows(rows)

    def test_full_green_capture_rejects_nonbijective_basis_map(self):
        rows=[]
        labels=[(0,0,0),(0,0,1),(1,-1,0),(1,0,0),(1,1,0),
                (1,-1,1),(1,0,1),(1,0,1)]
        for i,(ell,m,spin) in enumerate(labels):
            rows.append({'kind':'basis_map','stage':'RS_GREEN','site':'0','row':str(i),
                         'l':str(ell),'col':str(m),'spin':str(spin)})
        for i in range(8):
            for j in range(8):
                rows.append({'kind':'green_function','stage':'RS_GREEN','site':'0','k':'0',
                             'row':str(i),'col':str(j)})
        with self.assertRaisesRegex(ValueError,'duplicate \(l,m,spin\) basis label'):
            c3.validate_full_green_rows(rows)

    def test_c3_has_only_registered_scientific_gates_and_inclusive_bounds(self):
        data={'target':{'carrier_by_spin_ry':[1.0,2.0],'q2prime_by_spin_ry2':[0.0,0.2],
                        'capture_q2prime_by_spin_ry2':[0.0,0.2]}}
        cases={1:{'sumev':{1:10.0,2:20.0}},2:{'sumev':{1:9.0,2:17.98}}}
        got=c3.summarize_c3(data,cases)
        self.assertEqual(got['gates'],{
            'carrier_ratio_each_spin_within_1_percent':True,
            'q2prime_each_spin_nonnegative':True,
            'all_c3_gates_pass':True})
        self.assertAlmostEqual(got['carrier_ratio_by_spin']['2']['ratio'],1.01,places=14)

    def test_capture_integrity_binds_c2_state_and_all_trace_rows(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            deck=root/'input.nml'; seed=root/'Fe_out.nml'; binary=root/'rslmto.x'; source=root/'self.f90'; driver=root/'driver'
            hashes={p.name:write(p,p.name.encode()) for p in (seed,binary,source,driver)}
            hashes[deck.name]=write(deck,b"&self\nuse_kspace=.false.\n/\n&control\ncalctype='B'\nnsp=1\n"
                                       b"recur='chebyshev'\n/\n&hamiltonian\nhubbard_radial_moments='local'\n/")
            lock={'schema':c3.LOCK_SCHEMA,'status':c3.LOCK_STATUS,
                  'system':{'backend':'RS','calctype':'B','recurrence':'chebyshev','radial_moments':'local',
                            'nsp':1,'lmax':2,'u_ev':2.0,'j_ev':0.0},
                  'acceptance':{'c2_gate_status':'PASS','c1_rs_status':'ACCEPTED','convergence':{
                      'outer_diff':1e-8,'last20_d_occupation_by_spin':{'1':[1.0]*20,'2':[2.0]*20},
                      'last20_sphere_moment_muB':[3.0]*20}},
                  'accepted_endpoint':{'run_id':'accepted','accepted_iteration':20,
                       'deck_path':str(deck),'deck_sha256':hashes[deck.name],
                       'checkpoint_path':str(seed),'checkpoint_sha256':hashes[seed.name]},
                  'executable':{'binary_path':str(binary),'binary_sha256':hashes[binary.name],
                                'source_path':str(source),'source_sha256':hashes[source.name]}}
            lockp=root/'lock.json'; lockp.write_text(json.dumps(lock))
            binding=c3.validate_c2_lock(lock,lockp)
            cap=root/'capture'; cap.mkdir()
            rows=[]
            def row(kind,stage,site='0',typ='0',ell='-1',spin='0',k='-1',r='-1',col='-1',order='0'):
                rows.append([kind,'capture-run','1.0','1',stage,site,typ,ell,spin,k,'-1',r,col,order,'0','0','0','0','0'])
            row('basis_map','RS_GREEN',ell='0',spin='0',r='0',col='0')
            row('basis_map','RS_GREEN',ell='0',spin='1',r='1',col='0')
            for i in range(2):
                for j in range(2):
                    row('green_function','RS_GREEN',k='0',r=str(i),col=str(j))
            for kind in c3.REQUIRED_KINDS:
                if kind in {'basis_map','green_function','radial_replay_mesh','radial_replay_ec_seed','radial_replay_eigen_seed','m0_carrier_input','m0_carrier_local','m0_carrier_control'}: continue
                row(kind,'RADIAL_REPLAY')
            for i in range(495): row('radial_replay_mesh','RADIAL_REPLAY',spin='0',k='-1',r=str(i))
            for i in range(10): row('radial_replay_ec_seed','RADIAL_REPLAY',r=str(i))
            for i in range(6): row('radial_replay_eigen_seed','RADIAL_REPLAY',k=str(i))
            for sp in ('1','2'):
                row('m0_carrier_input','FROZEN_REPLAY',spin=sp)
                row('m0_carrier_local','FROZEN_REPLAY',spin=sp)
                row('m0_carrier_control','VOUT_CONTROL',spin=sp)
                row('m0_carrier_control','P_TRANSPOSE_CONTROL',spin=sp)
            fields=c3.rb3.CSV_FIELDS
            csvp=cap/'records.csv'
            with csvp.open('w',newline='') as f:
                w=csv.writer(f); w.writerow(fields); w.writerows(rows)
            manifest={'schema':c3.TRACE_SCHEMA,'backend':'RS','csv_sha256':c3.sha(csvp),
                      'executable_sha256':hashes[binary.name],'source_sha':hashes[source.name],
                      'input_sha256':{'accepted_checkpoint':hashes[seed.name]},
                      'run_id':'capture-run','iteration':1}
            mfp=cap/'manifest.json'; mfp.write_text(json.dumps(manifest))
            rows_out, info=c3.validate_capture(cap,binding,source,binary,driver)
            self.assertEqual(info['run_id'],'capture-run')
            self.assertTrue(info['green_shape']['all_green_matrices_complete'])

    def test_c3_does_not_floor_zero_carriers_or_hide_negative_q2(self):
        data={'target':{'carrier_by_spin_ry':[0.0,1.0],'q2prime_by_spin_ry2':[-1e-15,0.0],
                        'capture_q2prime_by_spin_ry2':[-1e-15,0.0]}}
        cases={1:{'sumev':{1:1.0,2:2.0}},2:{'sumev':{1:1.0,2:1.0}}}
        got=c3.summarize_c3(data,cases)
        self.assertIsNone(got['carrier_ratio_by_spin']['1']['ratio'])
        self.assertFalse(got['gates']['carrier_ratio_each_spin_within_1_percent'])
        self.assertFalse(got['gates']['q2prime_each_spin_nonnegative'])


if __name__ == '__main__':
    unittest.main()
