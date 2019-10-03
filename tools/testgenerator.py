import os
import subprocess
import plotslice as ps
import pylab as p
import numpy as np
import shutil

class TestGenerator():
    size_Mpc = 100.0
    powerlaw = None
    base_ncells = 32
    zoom_factor = [4]
    zoom_ncells = [32]
    modification_name = ['overdensity']
    modification_value = [0.15]
    modification_pos = [[size_Mpc/2]*3]
    modification_size = [None]
    path_to_IC = os.path.join(os.path.dirname(os.path.dirname(__file__)), "IC")

    setup_block = """
    Om  0.279
    Ol  0.721
    s8  0.817
    zin	99  
    seedfourier_parallel	8896131
    """

    @property
    def powspec_block(self):
        if self.powerlaw is None:
            return """camb	../camb_transfer_kmax40_z0.dat"""
        else:
            return """ns %.2f
            powerlaw_amplitude 1.0
            """%self.powerlaw

    output_block = """
    outname test
    outdir ./
    outformat 2"""

    @property
    def equiv_ncells(self):
        return self.zoom_ncells[-1] * np.prod(self.zoom_factor)

    @property
    def equiv_grid_block(self):
        return """basegrid %.2f %d
               zerolevel 0"""%(self.size_Mpc, self.equiv_ncells)

    @property
    def zoom_grid_block(self):
        s = "basegrid %.2f %d\n"%(self.size_Mpc, self.base_ncells)

        for fac, ncell in zip(self.zoom_factor, self.zoom_ncells):
            s+="""
                centre %.2f %.2f %.2f
                select_nearest
                zoomgrid %d %d
                """%(self.size_Mpc/2, self.size_Mpc/2, self.size_Mpc/2,
                     fac, ncell)

        for i in range(len(self.zoom_ncells)+1):
            s += "zerolevel %d\n" % i
        return s

    @property
    def modification_block(self):
        block=""
        for size, pos, name, val in zip(self.modification_size, self.modification_pos,
                                        self.modification_name, self.modification_value):
            if size:
                select_line = "select_sphere %.2f"%(size)
            else:
                select_line = "select_nearest"
            block+= """centre %.2f %.2f %.2f
            %s
            modify %s absolute %f
            """%(pos[0], pos[1], pos[2], select_line,
                                     name, val)
        return block

    @property
    def equiv_finalise_block(self):
        return """done
        dump_grid 0
        """

    @property
    def zoom_finalise_block(self):
        s = self.equiv_finalise_block
        for i in range(1,len(self.zoom_ncells)+1):
            s += "dump_grid %d\n" % i
        return s

    def dir_name(self, for_zoom):
        if for_zoom:
            return "zoom"
        else:
            return "equiv"

    def paramfile_name(self, for_zoom):
        return self.dir_name(for_zoom)+"/paramfile.txt"

    def _strip_whitespace(self, s: str):
        new_s = ""
        for l in s.splitlines():
            new_s+=l.strip()+"\n"
        return new_s

    def write_paramfile(self, for_zoom):
        grid_block = self.zoom_grid_block if for_zoom else self.equiv_grid_block
        finalise_block = self.zoom_finalise_block if for_zoom else self.equiv_finalise_block
        with open(self.paramfile_name(for_zoom), 'w') as f:
            print(self._strip_whitespace(self.setup_block), file=f)
            print(self._strip_whitespace(self.powspec_block), file=f)
            print(self._strip_whitespace(self.output_block), file=f)
            print(self._strip_whitespace(grid_block), file=f)
            print(self._strip_whitespace(self.modification_block), file=f)
            print(self._strip_whitespace(finalise_block), file=f)

    def setup(self, for_zoom):
        try:
            os.mkdir(self.dir_name(for_zoom))
        except FileExistsError:
            pass
        self.write_paramfile(for_zoom)

    def cleanup(self):
        for zoom in True, False:
            try:
                shutil.rmtree(self.dir_name(zoom))
            except FileNotFoundError:
                pass


    def run(self, for_zoom):
        subprocess.run([self.path_to_IC, self.paramfile_name(for_zoom)]).check_returncode()

    def run_all(self):
        for zoom in True, False:
            self.setup(zoom)
            self.run(zoom)

    def make_plots(self):
        p.set_cmap('RdBu_r')
        ax = p.subplot(221)
        p.title("Modification")
        for zoom in True, False:
            ps.plot1dslice(self.dir_name(zoom)+"/", slice_y=self.modification_pos[0][1], slice_z=self.modification_pos[0][2])
        new_ax = p.subplot(222, sharex=ax)
        ps.plot1dslice(self.dir_name(True)+"/", slice_y=self.modification_pos[0][1], slice_z=self.modification_pos[0][2],
                       diff_prefix=self.dir_name(False)+"/",vmin=-0.15,vmax=0.15)

        p.title("Relative error")
        p.ylim(-0.05,0.05)
        new_ax.yaxis.set_label_position('right')
        new_ax.yaxis.tick_right()

        ax = p.subplot(223)

        ps.plotslice(self.dir_name(True)+"/", slice=self.modification_pos[0][2])
        ax.text(0.05, 0.05, r"Colour scale: overdensity $\pm 0.15$",transform=ax.transAxes)
        new_ax = p.subplot(224, sharex=ax, sharey=ax)
        new_ax.yaxis.set_label_position('right')
        new_ax.yaxis.tick_right()
        ps.plotslice(self.dir_name(True)+"/", diff_prefix=self.dir_name(False)+"/", slice=self.modification_pos[0][2],
                     vmin=-0.05,vmax=0.05)
        new_ax.text(1.05, 0.05, r"$\pm 5\%$ of peak",transform=ax.transAxes)


    def go(self):
        self.cleanup()
        self.run_all()
        self.make_plots()

class VelTestGenerator(TestGenerator):
    modification_name = 'vx'
    modification_value = 5000