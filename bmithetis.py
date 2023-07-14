"""
Basic Model Interface (BMI) implementation for Thetis 2-D configuration
"""
from bmipy import Bmi
from colours import clrtxt
from firedrake import *

import numpy as np
from thetis import *
import time as time_mod
import warnings
warnings.simplefilter(action="ignore", category=DeprecationWarning)



class BmiThetis(Bmi):
    """
    Basic Model Interface implementation of 2-D Thetis setup
    """
    def __init__(self):
        """
        Create a BMI-refactored 2-D Thetis model that is ready for initialisa-
        tion
        """
        # Specify the name of the model component
        self._name = "2-D Thetis model"

        # 1-D array containing the names of the variables that the model can
        # use from other models implementing BMI, i.e. the names of the input
        # variables
        self._input_var_names = [
            "sea_surface_water_wave__significant_height",
            "sea_surface_water_wave__wavelength",
            "sea_surface_water_wave__direction",
            "sea_surface_water_wave__breaking_fraction"
        ]

        # 1-D array containing the names of the variables that the model can
        # provide to ther model implementing BMI, i.e. the names of the output
        # variables
        self._output_var_names = [
            "sea_water_surface__elevation",
            "sea_water_flow__x_component_of_velocity",
            "sea_water_flow__y_component_of_velocity",
            "sea_water_flow__x_and_y_component_of_velocity"
        ]

        # Dictionary containing the grid identifiers as keys and the correspon-
        # ding grid types as values
        self._grid_types = {
            0 : "uniform_rectilinear",
            1 : "rectilinear",
            2 : "structured_quadrilateral",
            3 : "unstructured"
        }

        # The data type of the variable as it's stored in memory by the model
        self._var_type = "float"

        # Dictionary containing the units of the given variable, for both input
        # and output variables
        self._var_units = {
            "sea_surface_water_wave__significant_height" : "m",
            "sea_surface_water_wave__wavelength" : "m",
            "sea_surface_water_wave__direction" : "degrees",
            "sea_surface_water_wave__breaking_fraction" : "None",
            "sea_water_surface__elevation" : "m",
            "sea_water_flow__x_component_of_velocity" : "m s-1",
            "sea_water_flow__y_component_of_velocity" : "m s-1",
            "sea_water_flow__x_and_y_component_of_velocity" : "m s-1"

        }

        # Time unit as employed by Thetis
        self._time_units = "s"

        # Not sure what this is
        self._thetis_names = {
            "sea_water_surface__elevation" : "elev",
            "sea_water_flow__x_component_of_velocity": "uv",
            "sea_water_flow__y_component_of_velocity": "uv",
            "sea_water_flow__x_and_y_component_of_velocity": "uv"
        }


    def initialize(self, filename=None):
        """
        Perform all tasks that are to take place before entering the model's
        time loop. These are:
        * Read the configuration file
        * Configure the mesh
        * Get the mesh coordinates
        * Perform any coordinates conversions
        * Define the necessary function spaces
        * Configure the bathymetry
        * Define the viscosity field
        * Define the manning coefficient field
        * Define the coriolis parameter field
        * Define a TPXO tide interpolator
        * Define any functions requested by the user
        * Setup the Shallow-Water Equation Boundary Conditions
        * Get the MPI status
        * Setup the Wave-Current Interactions
        * Create a flowsolver
        Inputs:
        self     :
        filename : The filepath of the json for the simulation configuration

        Outputs:
        None

        Returns:
        None
        """
        from bmi_tools import BathymetryTools, SetupTools, ThetisTools
        from bmi_tools import WCITools, CoordsTools, TideTools, MPITools
        from bmi_tools import IOTools, WindTools, BCTools, FinaliseTools
        from bmi_tools import ForcingTools, WEoCTools
        from pathlib import Path

        # Load configuration file
        IOTools.load_json_file(self, filename)

        # Set coupling and and simulation status
        SetupTools.set_status(self)

        # Load mesh filepath
        SetupTools.configure_mesh(self)

        # Get mesh coordinates
        CoordsTools.get_mesh_coordinates(self)

        # Check if there is any conversion of coordinates that needs to be done
        if "coordinates conversion" in self.config:
            CoordsTools.convert_coordinates(self)

        # Check if the simulation is run in parallel
        MPITools.get_mpi_status(self)

        # Define function spaces
        SetupTools.define_function_spaces(self)

        # Setup bathymetry
        BathymetryTools.bathymetry_field(self)

        # Setup eddy viscosity
        if "viscosity parameters" in self.config:
            SetupTools.viscosity_field(self)

        # Setup manning coefficient
        if "manning parameters" in self.config:
            SetupTools.manning_field(self)

        # Account for Coriolis forcing
        if "coriolis parameters" in self.config:
            SetupTools.coriolis_field(self)

        # Functions (for boundary conditions and other)
        if "function definition" in self.config:
            SetupTools.user_defined_function(self)



        ## TIde forcings
        if "tide forcing" in self.config:
            TideTools.configure_tide_forcing(self)
        # if "TPXO forcing" in self.config:
        #     TideTools.configure_tide_interp(self)

        ## Velocity forcing
        if "uv forcing" in self.config:
            ForcingTools.configure_uv_interp(self)

        if "wind forcing" in self.config:
            WindTools.wind_field(self)

        # Boundary conditions
        SetupTools.swe_boundary_conditions(self)



        # Wave-Current Interactions (WCI)
        if "WCI" in self.config:
            WCITools.configure_wci(self)

        # Wave-Effects on Currents (Offline Coupling)
        if "WEoC" in self.config:
            WEoCTools.configure_weoc_forcing(self)

        # Configure 2-D flowsolver for Shallow-Water Equations
        ThetisTools.create_flowsolver(self)

        # Check the update forcings
        BCTools.check_user_defined_forcings(self)

        # Check theparameters for finalising the model
        FinaliseTools.check_finalise_params(self)

        # If the coupling is both ways or from SWAN to Thetis
        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "Thetis-to-SWAN"
        ):
            WCITools.calculate_output_fields(self)

        # Get the timestep
        self.dt = self.get_time_step()



        return


    def update(self):
        """
        Advance the model by a single timestep. This means that we update the
        Wave-Effects-On-Currents objects if we account for WCI (Wave-Current
        Interactions), advance the solution and update the time and the relev-
        ant indices, as well as calculate any callbacks, and export any output.
        Inputs:
        self :

        Outputs:
        None

        Returns:
        None
        """
        from bmi_tools import ThetisTools
        def update_forcings(t):
            from bmi_tools import BCTools
            # Check that the user has specified boundary conditions to be upda-
            # ted
            if "update forcings" in self.config:
                BCTools.update_user_defined_forcings(self, t)
            else:
                message = clrtxt.bold + clrtxt.darkcyan + "Attention: " + \
                    clrtxt.end + "No forcing are updated"
                raise SystemError(message)


        # Calculate ramp coefficient
        ThetisTools.calculate_ramp_coeff(self, self.solver_obj.simulation_time)

        # Get the current time is seconds
        cputimestamp = time_mod.perf_counter()

        # If the coupling is both ways or from SWAN to Thetis
        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "SWAN-to-Thetis" or
                "WEoC" in self.config
        ):
            # Update the Wave Effects On Currents
            self.weoc.update(
                hs_wave=self.hs,
                dir_wave=self.dir,
                l_wave=self.wlen,
                qb=self.qb,
                c1=self.ramp_coeff
            )

        # Advance the equation for one timestep
        self.solver_obj.timestepper.advance(
            self.solver_obj.simulation_time,
            update_forcings=update_forcings
        )

        ## Move to next timestep
        # Update the iteration count
        self.solver_obj.iteration += 1
        # Update the simulation time
        self.solver_obj.simulation_time += self.solver_obj.dt
        # EValuate any callbacks
        self.solver_obj.callbacks.evaluate(mode='timestep')

        # Write the solution to the stdout file
        if self.solver_obj.iteration in self.export_ind:
            # Update the export index
            self.solver_obj.i_export += 1
            # Calculate how much time advancing the model to the next timestep
            # took
            cputime = time_mod.perf_counter() - cputimestamp
            # Get the current time in seconds
            cputimestamp = time_mod.perf_counter()
            # Print a summary of the model state on the stdout
            self.solver_obj.print_state(cputime)
            # Export all fields to disk, i.e. hdf5, pvd, etc
            self.solver_obj.export()

        return


    def update_until(self, t_end):
        """
        Update Thetis to a particular time. Depending on the time specified,
        the method update() only may be called or in collaboration with Thetis
        inherent method solver_obj.iterate_coupled. If we account for Wave-Cur-
        rent Interactions, we will also calculate the values for the export fi-
        elds.
        Inputs:
        self :
        t_end : A float containing the time
        """
        from bmi_tools import WCITools

        def update_forcings(t):
            from bmi_tools import BCTools, ThetisTools
            ThetisTools.calculate_ramp_coeff(self, t)

            # Check that the user has specified boundary conditions to be upda-
            # ted
            if "update forcings" in self.config:
                BCTools.update_user_defined_forcings(self, t)
            else:
                message = clrtxt.bold + clrtxt.darkcyan + "Attention: " + \
                          clrtxt.end + "No forcing are updated"
                raise SystemError(message)

            ## Update Wave Current INteractions
            if (
                    self.coupl_stat == "2-way" or
                    self.coupl_stat == "SWAN-to-Thetis"
            ):
                # If the wave fields/characteristics needs to be interpolated
                if self.wave_char_interp:
                    ## Calculate the new value of wave characteristics
                    # For the significant wave height
                    self.hs.interpolate(
                        self.hs_diff * Constant(self.wave_interp_counter) + \
                        self.hs_old
                    )
                    # Fort the mean wave direction
                    self.dir.interpolate(
                        self.dir_diff * Constant(self.wave_interp_counter) + \
                        self.dir_old
                    )
                    # For the percentage of wave breaking
                    self.qb.interpolate(
                        self.qb_diff * Constant(self.wave_interp_counter) + \
                        self.qb_old
                    )
                    # For the mean wavelength
                    self.wlen.interpolate(
                        self.wlen_diff * Constant(self.wave_interp_counter) + \
                        self.wlen_old
                    )

                    # Update counter of wave characteristic interpolation
                    self.wave_interp_counter = self.wave_interp_counter + 1

                    # Update the Wave Effects On Currents
                    self.weoc.update(
                        hs_wave=self.hs,
                        dir_wave=self.dir,
                        l_wave=self.wlen,
                        qb=self.qb,
                        c1=self.ramp_coeff
                    )

            if "WEoC" in self.config:
                # Update the Wave Effects On Currents
                self.weoc.update(
                    hs_wave=self.hs,
                    dir_wave=self.dir,
                    l_wave=self.wlen,
                    qb=self.qb,
                    c1=self.ramp_coeff
                )

        ## Specify the wave characteristics for the 1st update
        # If the coupling is both ways or from SWAN to Thetis
        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "SWAN-to-Thetis"
        ):
                # If the wave fields need to be interpolates
            if self.wave_char_interp:
                # STart counter
                self.wave_interp_counter = 1
                # Calculate the "step" of each wave characteristic at each
                # point
                self.hs_diff.interpolate(
                    Constant(1 / self.no_wave_interp) * \
                    (self.hs_new - self.hs_old)
                )
                self.dir_diff.interpolate(
                    Constant(1 / self.no_wave_interp) * \
                    (self.dir_new - self.dir_old)
                )
                self.wlen_diff.interpolate(
                    Constant(1 / self.no_wave_interp) * \
                    (self.wlen_new - self.wlen_old)
                )
                self.qb_diff.interpolate(
                    Constant(1 / self.no_wave_interp) * \
                    (self.qb_new - self.qb_old)
                )

                # Calculate the new value of wave characteristics
                self.hs.interpolate(self.hs_old + self.hs_diff)
                self.dir.interpolate(self.dir_old + self.dir_diff)
                self.qb.interpolate(self.qb_old + self.qb_diff)
                self.wlen.interpolate(self.wlen_old + self.wlen_diff)

                # Update counter
                self.wave_interp_counter = self.wave_interp_counter + 1

            # If Thetis dt = coupling dt
            else:
                self.hs = self.hs_new
                self.dir = self.dir_new
                self.wlen = self.wlen_new
                self.qb = self.qb_new

        self.update()

        ## Check if the current simulation time
        # Calculate the time difference between the current simulation time and
        # the time provided as input t_end
        t_diff = abs(self.solver_obj.simulation_time - t_end)
        # If the time difference isn't negligible
        if t_diff > 10**(-5):
            self.solver_obj.options.simulation_end_time = float(t_end)
            self.solver_obj.iterate_coupled(update_forcings=update_forcings)

        # If the coupling is both ways or from SWAN to Thetis
        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "Thetis-to-SWAN"
        ):
            WCITools.calculate_output_fields(self)

        # If SWAN passes information to Thetis
        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "SWAN-to-Thetis"
        ):
            # Pass the "new" value, i.e. the final value to the old array
            self.hs_old.interpolate(self.hs_new)
            self.dir_old.interpolate(self.dir_new)
            self.qb_old.interpolate(self.qb_new)
            self.wlen_old.interpolate(self.wlen_new)

        self.comm.Barrier()

        return


    def finalize(self):
        """
        Add description
        """
        from bmi_tools import FinaliseTools, IOTools
        from pathlib import Path

        print_output("- Exporting final state of Thetis")

        # Check if the user has defined any time information to be included in
        # the file names
        time_info = FinaliseTools.check_time_info(self)

        # Load the output directory
        output_dir = self.config["finalize parameters"]["output directory"]
        Path(output_dir).mkdir(exist_ok=True, parents=True)

        if time_info is None:
            IOTools.export_h5(
                self, self.uv, "velocity", output_dir + "velocity"
            )
            IOTools.export_h5(
                self, self.elev, "elevation", output_dir + "elevation"
            )
            # If Thetis accepts information from SWAN
            if (
                    (self.coupl_stat == '2-way') or
                    (self.coupl_stat == "SWAN-to-Thetis")
            ):
                # Export the significant waveheight
                IOTools.export_h5(
                    self,
                    self.hs,
                    "hs",
                    output_dir + "hs"
                )
                # Export the mean wave direction
                IOTools.export_h5(
                    self,
                    self.dir,
                    "dir",
                    output_dir + "dir"
                )
                # Export the mean wavelength
                IOTools.export_h5(
                    self,
                    self.wlen,
                    "wlen",
                    output_dir + "wlen"
                )
                # Ecport the percentage of wave breaking
                IOTools.export_h5(
                    self,
                    self.qb,
                    "qb",
                    output_dir + "qb"
                )

        else:
            IOTools.export_h5(
                self, self.uv, "velocity", output_dir + "velocity_" + time_info
            )
            IOTools.export_h5(
                self,
                self.elev,
                "elevation",
                output_dir + "elevation_" + time_info
            )
            # If Thetis accepts information from SWAN
            if (
                    (self.coupl_stat== '2-way') or
                    ( self.coupl_stat == "SWAN-to-Thetis")
            ):
                # Export the significant waveheight
                IOTools.export_h5(
                    self,
                    self.hs,
                    "hs",
                    output_dir + "hs_" + time_info
                )
                # Export the mean wave direction
                IOTools.export_h5(
                    self,
                    self.dir,
                    "dir",
                    output_dir + "dir_" + time_info
                )
                # Export the mean wavelength
                IOTools.export_h5(
                    self,
                    self.wlen,
                    "wlen",
                    output_dir + "wlen_" + time_info
                )
                # Ecport the percentage of wave breaking
                IOTools.export_h5(
                    self,
                    self.qb,
                    "qb",
                    output_dir + "qb_" + time_info
                )




    def get_component_name(self):
        """
        Get the name of the model component as a string

        Outputs:
        _name : Name of model component
        """
        return self._name


    def get_input_item_count(self):
        """
        The number of variables the model can use from other models

        Outputs:
        len(self._input_var_names) : Maximum number of imported variables
        """
        return len(self._input_var_names)


    def get_output_item_count(self):
        """
        The number of variables the model can provide to other models

        Outputs:
        len(self._output_var_names) : Maximum number of exported variables
        """
        return len(self._output_var_names)


    def get_input_var_names(self):
        """
        Names of the variables model can use from other models

        Outputs:
        self._input_var_names : List of imported variables names
        """
        return self._input_var_names


    def get_output_var_names(self):
        """
        Names of the variables the model can provide to other models

        Outputs:
        self._output_var_names : List of exported variables names
        """
        return self._output_var_names


    def get_var_grid(self):
        """

        """

    def get_var_type(self, var_name):
        """
        Get the data type of one of the model's variables

        Inputs:
        var_name : The name of the variable.

        Output:
        self._var_type : The type of the var_name variable
        """
        # Check that the variable name is valid
        if (
            (var_name in self._input_var_names) or
            (var_name in self._output_var_names)
        ):
            return self._var_type
        else:
            print_output(f"Variable name '{var_name}' is not valid. "\
                "Please select one of the following:")
            print_output("- For model's input variables:")
            print_output("\n".join(self._input_var_names))
            print_output("- For model's output variables:")
            print_output("\n".join(self._output_var_names))
            exit()


    def get_var_units(self, var_name):
        """
        Get the units of one of the model's variables.

        Inputs:
        var_name : The name of the variable.

        Output:
        unit : The unit of the var_name variable
        """
        # Check that the variable name is valid
        if (
            (var_name in self._input_var_names) or
            (var_name in self._output_var_names)
        ):
            unit = self._var_units[var_name]
            return unit
        else:
            print_output(f"Variable name '{var_name}' is not valid. "\
                "Please select one of the following:")
            print_output("- For model's input variables:")
            print_output("\n".join(self._input_var_names))
            print_output("- For model's output variables:")
            print_output("\n".join(self._output_var_names))
            exit()


    def get_var_itemsize(self):
        """
        """

    def get_var_nbytes(self):
        """
        """

    def get_var_location(self):
        """
        """

    def get_current_time(self):
        """
        Get the current time of the model

        Outputs:
        self.solver_obj.simulation_time : Current time of the model [s]
        """
        return self.solver_obj.simulation_time


    def get_start_time(self):
        """
        Get the start time of the model
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.t_start : Start time of the model [s]
        """
        return self.t_start


    def get_end_time(self):
        """
        Get the end time of the model
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.t_end : End time of the model [s]
        """
        return self.t_end


    def get_time_units(self):
        """
        Get the time units used in the model

        Outputs:
        self._time_units : Time units used in the model [str]
        """
        return self._time_units


    def get_time_step(self):
        """
        Get the time step used in the Thetis model

        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.solver_obj.options.timestep : Model's time-step [s]
        """
        return self.solver_obj.options.timestep


    def get_value(self, var_name):
        """
        Get variable values

        Inputs:
        var_name : Variable name [str]
        grid_name : The grid of the variable
        """
        from bmi_tools import WCITools

        if var_name not in self._output_var_names:
            raise KeyError(var_name)

        src = WCITools.get_output_field(self, var_name)
        return src

    def get_value_ptr(self):
        """
        """

    def get_value_at_indices(self):
        """
        """

    def set_value(self, var_name, new_vals):
        """

        """
        from bmi_tools import WCITools

        if var_name not in self._input_var_names:
            raise KeyError(var_name)

        WCITools.set_input_field(self, var_name, new_vals)

    def set_value_at_indices(self):
        """
        """

    def get_grid_rank(self):
        """
        """

    def get_grid_size(self):
        """
        """

    def get_grid_type(self):
        """
        """

    def get_grid_shape(self):
        """
        """

    def get_grid_spacing(self):
        """
        """

    def get_grid_origin(self):
        """
        """

    def get_grid_x(self):
        """
        """

    def get_grid_y(self):
        """
        """

    def get_grid_z(self):
        """
        """

    def get_grid_node_count(self):
        """
        """

    def get_grid_edge_count(self):
        """
        """

    def get_grid_face_count(self):
        """
        """

    def get_grid_edge_nodes(self):
        """
        """

    def get_grid_face_edges(self):
        """
        """

    def get_grid_face_nodes(self):
        """
        """

    def get_grid_nodes_per_face(self):
        """
        """
