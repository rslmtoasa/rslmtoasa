"""Rotate the d-shell spin frame of a collinear (z) surface seed to a new axis: spinor P' = W P W^+, W = exp(-i th/2 n.sigma) x 1_5,
spinor_d = 1 blocks written into <label>_dm.nml and mom updated in the potential files.
Usage: rot_seed.py <dir> <theta_deg> <label...>   (rotation about y by theta; z -> x for +90)"""
import sys,re,numpy as np
d=sys.argv[1]; th=np.deg2rad(float(sys.argv[2])); labels=sys.argv[3:]
sy=np.array([[0,-1j],[1j,0]]); W=np.kron(np.cos(th/2)*np.eye(2)-1j*np.sin(th/2)*sy,np.eye(5))
n=[np.sin(th),0.0,np.cos(th)]
for l in labels:
    p=f"{d}/{l}_dm.nml"; txt=open(p).read(); P=np.zeros((10,10),complex)
    for s in (1,2):
        for t in (1,2):
            for i in range(1,6):
                a=re.search(rf"dm_re\(3,{s},{t},{i},:\) =(.*)",txt); b=re.search(rf"dm_im\(3,{s},{t},{i},:\) =(.*)",txt)
                if a: P[(s-1)*5+i-1,(t-1)*5:t*5]=np.array(a[1].split(),float)+1j*np.array(b[1].split(),float)
    Pn=W@P@W.conj().T
    txt=re.sub(r"(?m)^ dm_(re|im)\(3,.*\n","",txt)
    new=" spinor_d = 1\n"
    for s in (1,2):
        for t in (1,2):
            for i in range(1,6):
                blk=Pn[(s-1)*5+i-1,(t-1)*5:t*5]
                new+=f" dm_re(3,{s},{t},{i},:) = "+"  ".join(f"{v:.16E}" for v in blk.real)+"\n"
                new+=f" dm_im(3,{s},{t},{i},:) = "+"  ".join(f"{v:.16E}" for v in blk.imag)+"\n"
    txt=re.sub(r"(?m)^ spinor_d = 1\n","",txt); txt=txt.replace("\n/","\n"+new+"/",1) if txt.rstrip().endswith("/") else txt+new
    open(p,"w").write(txt)
    q=f"{d}/{l}.nml"; t2=open(q).read(); t2=re.sub(r"(?m)^ hubbard_nsp = .*$"," hubbard_nsp = 2",t2); open(q,"w").write(t2)
for l in ("ES",)+tuple(labels):
    q=f"{d}/{l}.nml"; t2=open(q).read(); t2=re.sub(r"(?m)^ mom = .*$"," mom = %.16f, %.16f, %.16f"%tuple(n),t2); t2=re.sub(r"(?m)^ hubbard_axis = .*$"," hubbard_axis = %.16f, %.16f, %.16f"%tuple(n),t2); open(q,"w").write(t2)
print("rotated",labels,"n=",n)
