#!/usr/bin/env python3
"""Source-grounded namelist and recipe check for the proposed Ni U0 bulk host."""
import hashlib, json, sys
from pathlib import Path
import f90nml
sys.path.insert(0,'/private/tmp/c52-c-private')
from fortran_namelist_schema import source_schema, validate_input
REPO=Path('/Users/andersb/Jobb/rslmto_fable_v4b')
ROOT=Path('/private/tmp/c52-rd-ni-u0-host')
SCHEMA=source_schema(REPO)
paths=[ROOT/'immutable/input.nml',ROOT/'immutable/Ni-imp.nml',ROOT/'immutable/Ni-imp_out.nml']
validations=[validate_input(p,SCHEMA) for p in paths]
n=f90nml.read(str(paths[0]))
seed=f90nml.read(str(paths[2]))
def req(ok,msg):
 if not ok: raise SystemExit('preflight assertion failed: '+msg)
req(n['calculation']['pre_processing']=='bravais','calculation bravais')
lat=n['lattice']; req(lat['rc']==80 and lat['crystal_sym']=='fcc','fcc rc80')
req(abs(lat['alat']-3.52340)<1e-10 and abs(lat['wav']-1.37692)<1e-10,'accepted Ni lattice constants')
req(lat['ntype']==1 and abs(lat['ct'][0]-4.0)<1e-12 and abs(lat['r2']-16.0)<1e-12,'single Ni bulk site cluster settings')
req(n['atoms']['label'][0]=='Ni-imp' and n['atoms']['database']=='./','loader alias')
req(n['self']['use_kspace'] is False,'RS self.use_kspace=false')
req(n['control']['calctype']=='B' and n['control']['nsp']==1 and n['control']['lld']==21 and n['control']['recur']=='block' and n['control']['fresh_start'] is False,'bulk scalar block recipe and reload')
req(n['self']['nstep']==100 and abs(n['self']['conv_thr']-1e-30)<1e-45,'initial bounded capture cap')
e=n['energy']; req(e['fix_fermi'] is True and abs(e['fermi']+0.061141)<1e-12 and abs(e['energy_min']+1.0)<1e-12 and abs(e['energy_max']-0.5)<1e-12 and e['channels_ldos']==2500,'matched fixed-EF and energy window')
m=n['mix']; req(abs(m['beta']-0.1)<1e-12 and m['mixtype']=='broyden','matched Broyden mixing')
h=n['hamiltonian']; req(h['hubbard_radial_moments']=='local','explicit local radial-moment mode')
# f90nml may represent the sparse 2D Fortran arrays with only initialized bounds.
def cell(x,i,j):
 try: return float(x[i][j])
 except (IndexError,TypeError): return None
req(float(h['hubbard_u_general'][0][0])==0.0 and float(h['hubbard_j_general'][0][0])==0.0,'explicit sparse U=0/J=0 assignments (1,3) parsed')
req(seed['element']['symbol']=='Ni-imp' and seed['par']['lmax']==2,'accepted impurity seed type and lmax')
req(all(float(v)==0.0 for v in seed['par']['hubbard_u']) and all(float(v)==0.0 for v in seed['par']['hubbard_j']),'accepted seed stores U=J=0')
req(hashlib.sha256(paths[2].read_bytes()).hexdigest()=='addd6667159c9e4a5b6e4d3ff575143281877fb0b3c78832dccead1a34735e2d','accepted seed SHA')
out={'status':'PASS_SOURCE_SCHEMA_AND_RECIPE_ASSERTIONS','python':sys.executable,'f90nml':getattr(f90nml,'__version__','unknown'),'schema_source_root':str(REPO),'schema_source_files':[str(p.relative_to(REPO)) for p in sorted((REPO/'source/include_codes/namelists').glob('*.f90'))],'validated_decks':validations,'recipe':{'calctype':'B','system':'fcc Ni', 'rc':80,'lld':21,'recur':'block','nsp':1,'U_eV':0.0,'J_eV':0.0,'hubbard_radial_moments':'local','explicit_sparse_assignments':['hubbard_u_general(1,3)=0 eV','hubbard_j_general(1,3)=0 eV'],'EF_Ry':-0.061141,'fix_fermi':True,'energy_window_Ry':[-1.0,0.5],'beta':0.1,'mixtype':'broyden','initial_nstep':100,'conv_thr_Ry':1e-30,'reload_alias':'Ni-imp_out.nml','fresh_start':False},'seed_sha256':hashlib.sha256(paths[2].read_bytes()).hexdigest()}
print(json.dumps(out,indent=2))
