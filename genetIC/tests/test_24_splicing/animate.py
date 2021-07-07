import subprocess
import plotslice as ps
import pylab as p

def run_genetic(sphere_size):
    st = """
Om  0.279
Ol  0.721
#Ob  0.04
s8  0.817
zin	99

random_seed_real_space	8896131
camb	../camb_transfer_kmax40_z0.dat

outname test_1
outdir	 ./
outformat tipsy


basegrid 50.0 128

center 25 25 25
select_sphere %.1f
splice 8896132

done

dump_ps 0
dump_grid 0"""%sphere_size

    with open("paramfile_anim.txt","w") as f:
        f.write(st)

    subprocess.run(["../../genetIC", "paramfile_anim.txt"])

for frame in range(100):
    size = frame/5
    run_genetic(size)
    p.clf()
    ps.plotslice("./")
    p.savefig("anim/%.3d.png"%frame)
    ps.plotslice("./",diff_prefix="comp/",diff_normalized=False)
    p.savefig("anim_diff/%.3d.png"%frame)
