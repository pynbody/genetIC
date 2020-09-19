import os
import subprocess
import plotslice as ps
import pylab as p
import numpy as np
import shutil

_genetIC_folder = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "genetIC")

class TestGenerator():
    """Flexible test generator class.

    Example usage:

     x = TestGenerator()
     x.go()

    This runs genetIC twice, once with actual zoom regions and another time with an equivalent uniform-res grid
    then plots a four-panel plot to show the result and differences.

    To customise, you can either derive a class and change its settings or manually manipulate the object, e.g.
    by setting x.zoom_factor = [2,2]; x.zoom_ncells = [64, 64] you will get a double zoom.

    For all parameters that can be varied look at the top of the TestGenerator class.
    """

    # Here are the parameters you may wish to vary:
    test_name = "default" # used for plotting
    size_Mpc = 100.0 # size of the base box
    powerlaw = None # if None, use the test CAMB power spectrum; otherwise specifies a powerlaw power spectrum
    base_ncells = 64 # number of cells on the base grid
    zoom_factor = [4] # list of zoom factors (one for each level)
    zoom_ncells = [64] # list of zoom grid resolutions (one for each level)
    modification_name = ['overdensity'] # list of modification types (one for each modification)
    modification_value = [0.15] # list of modification target values (one for each modification)
    modification_pos = [[size_Mpc/2]*3] # list of modification target positions (one for each modification)
    zoom_cen = [size_Mpc/2, size_Mpc/2, size_Mpc/2] # centre of the zoom region
    modification_size = [None] # list of modification region sphere sizes (one for each modif, or None for point)
    path_to_IC =  os.path.join(_genetIC_folder,"genetIC") # path to the genetIC executable

    setup_block = """
    Om  0.279
    Ol  0.721
    s8  0.817
    zin	99
    random_seed	8896131
    """

    @property
    def powspec_block(self):
        if self.powerlaw is None:
            return """camb	"""+os.path.join(_genetIC_folder, "tests", "camb_transfer_kmax40_z0.dat")
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
                """%(self.zoom_cen[0],self.zoom_cen[1],self.zoom_cen[2],
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
        dump_vx 0
        """

    @property
    def zoom_finalise_block(self):
        s = self.equiv_finalise_block
        for i in range(1,len(self.zoom_ncells)+1):
            s += "dump_grid %d\n" % i
            s += "dump_vx %d\n" % i
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
        """Create the paramfile for either the zoom run or the equivalent uniform-res run"""
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
        """Create a folder and paramfile for either the zoom run or the equivalent uniform-res run"""
        try:
            os.mkdir(self.dir_name(for_zoom))
        except FileExistsError:
            pass
        self.write_paramfile(for_zoom)

    def cleanup(self):
        """Remove files created during a test run"""
        for zoom in True, False:
            try:
                shutil.rmtree(self.dir_name(zoom))
            except FileNotFoundError:
                pass


    def run(self, for_zoom):
        """Run either the zoom or the equivalent uniform-res"""
        subprocess.run([self.path_to_IC, self.paramfile_name(for_zoom)]).check_returncode()

    def run_all(self):
        """Set up and run both the zoom and equivalent ICs"""
        for zoom in True, False:
            self.setup(zoom)
            self.run(zoom)

    def make_plots(self, save_filename=None, vx=False):
        """Make a four-panel plot with 1d and 2d results and differences to the ideal case"""
        if vx:

            name = 'vx'
        else:

            name = 'grid'

        vmax = np.round(ps.get_peak_value(self.dir_name(False), name),2)
        vmin = -vmax

        if save_filename:
            p.figure(figsize=(9,9))
        if vx:
            p.set_cmap('RdBu_r')
        else:
            p.set_cmap('PuOr')
        ax = p.subplot(221)
        p.title("Modification")
        for zoom in True, False:
            ps.plot1dslice(self.dir_name(zoom)+"/", slice_y=self.modification_pos[0][1], slice_z=self.modification_pos[0][2],
                           name=name)
        new_ax = p.subplot(222, sharex=ax)
        ps.plot1dslice(self.dir_name(True)+"/", slice_y=self.modification_pos[0][1], slice_z=self.modification_pos[0][2],
                       diff_prefix=self.dir_name(False)+"/",name=name)

        p.title("Relative error")
        p.ylim(-0.05,0.05)
        new_ax.yaxis.set_label_position('right')
        new_ax.yaxis.tick_right()

        ax = p.subplot(223, sharex=ax)

        ps.plotslice(self.dir_name(True)+"/", slice=self.modification_pos[0][2],name=name,vmin=vmin,vmax=vmax)
        friendly_name = "overden" if name=="grid" else name
        ax.text(0.05, 0.05, r"Colour scale: %s $\pm %.2f$"%(friendly_name,vmax),transform=ax.transAxes)
        new_ax = p.subplot(224, sharex=ax, sharey=ax)
        new_ax.yaxis.set_label_position('right')
        new_ax.yaxis.tick_right()
        ps.plotslice(self.dir_name(True)+"/", diff_prefix=self.dir_name(False)+"/", slice=self.modification_pos[0][2],
                     vmin=-0.02,vmax=0.02,name=name)
        new_ax.text(1.05, 0.05, r"$\pm 2\%$ of peak",transform=ax.transAxes)

        if save_filename:
            p.subplots_adjust(wspace=0,hspace=0,left=0.1,right=0.9)
            p.savefig(save_filename,bbox_inches='tight')
            p.close()


    def go(self, save_filename=None):
        """Run genetIC for the zoom and equivalent uniform-res runs, then plot the results"""
        self.cleanup()
        self.run_all()
        self.make_plots(save_filename)

    def go_no_interaction(self):
        self.go("figures/"+self.test_name+"-overden-plot.pdf")
        self.make_plots("figures/"+self.test_name+"-vx-plot.pdf",True)

class VelTestGenerator(TestGenerator):
    test_name = 'velocity'
    modification_name = ['vx']
    modification_value = [5000]

class OffCentreTestGenerator(TestGenerator):
    test_name = 'offcen'
    modification_pos = [[TestGenerator.size_Mpc/2 + 3.0,
                         TestGenerator.size_Mpc/2, TestGenerator.size_Mpc/2]]

class OffCentreVelTestGenerator(OffCentreTestGenerator, VelTestGenerator):
    test_name = 'velocity-offcen'

class ZoomIsOffCentreVelTestGenerator(OffCentreTestGenerator, VelTestGenerator):
    test_name = 'velocity-grid-offcen'
    zoom_cen = [TestGenerator.size_Mpc/2 + 7.0, TestGenerator.size_Mpc/2, TestGenerator.size_Mpc/2]


class DoubleZoomTestGenerator(OffCentreTestGenerator):
    test_name = 'doublezoom'
    base_ncells = 64
    zoom_factor = [2, 2]
    zoom_ncells = [64, 64]

class DoubleModifTestGenerator(TestGenerator):
    modification_pos = TestGenerator.modification_pos + [[TestGenerator.size_Mpc/2 + 3.0,
                        TestGenerator.size_Mpc/2, TestGenerator.size_Mpc/2]]
    modification_size = TestGenerator.modification_size + [3.0]
    modification_value = TestGenerator.modification_value + [0.0]
    modification_name = ['overdensity']*2
    test_name = 'double_modif'

def _all_subclasses(cls):
    """Get all subclasses of cls (recursively)"""
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses__() for s in _all_subclasses(c)])

def all_runs():
    """Generate and run all known tests, saving output to a figures subfolder"""
    try:
        os.mkdir("figures")
    except:
        pass

    for C in _all_subclasses(TestGenerator):
        C().go_no_interaction()

if __name__ == "__main__":
    all_runs()
