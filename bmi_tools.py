"""
This file contains the following classes:
* CoordsTools
* BathymetryTools
* IOTools
* SetupTools
* ThetisTools
* DetectorTools
* WCITools
* TideTools
* MPITools
* BCTools
"""
from bmithetis import BmiThetis
from colours import clrtxt, ct
from firedrake import *
import numpy as np
import pandas as pd
from thetis import *


class CoordsTools(BmiThetis):
    """
    Mainly coordinate conversion in the case our mesh has different coordinates
    than the input files or it is necessary to procure results in another coor-
    dinate system
    Methods included in this class:
    * convert_coordinates(self)
    * convert_coords_to_another_epsg(self, dict, input_epsg, output_epsg)
    * convert_local_crs_to_latlon(self, dict, input_epsg, output_epsg)
    * get_mesh_coordinates(self)
    * get_mesh_coords_in_TPXO_crs(self)
    * get_mesh_coords_in_wind_crs(self, dict)
    * process_coords_conversion(self, dict)
    """

    def convert_coordinates(self):
        """
        Convert coordinates from one Coordinate Reference System (CRS) to ano-
        ther. Run some checks to verify the validity/accuracy of the options
        under the 'coordinates' entry in the json file. Once these checks have
        been completed, transform the coordinates
        Inputs:
        self

        Outputs:
        None

        Returns:
        self.dict_coords : A dictionary containing the name of the available
                           CRSes as keys. Each key corresponds to a dictionary
                           containing the following keys:
                           * x_coords : 1-D array containing the x-coordinates
                           * y_coords : 1-D array containing the y-coordinates
                           * x-coord name : A string containing the name of the
                                x-coordinates and their units, for example
                                'Longitude [o]', 'Easting [m]', etc
                           * y-coord name : A String containing the name of the
                                y-coordinates and their units, for example
                                'Latitude [o]', 'Northing [m], etc
        """
        # Get the dictionary specifying the coordinate conversion parameters
        dict = self.config["coordinates conversion"]
        # Get the number of Coordinate Reference System
        no_crs = int(dict["number of CRS"])
        # Get the number of coordinate conversion requested
        no_conv = int(np.shape(list(dict.keys()))[0] - 1)
        # Get the coordinate conversion requests
        conv = list(dict.keys())[1:]

        # Check that the number of CRS is more than 1
        if no_crs < 2:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "Although coordinate conversion has been " + \
                      "requested, the number of CRS to be available is " + \
                      clrtxt.bold + clrtxt.darkcyan + str(
                no_crs) + clrtxt.end + \
                      " instead of more than 1 as it should be. Please make sure" + \
                      " that either no coordinate conversion is required, i.e. " + \
                      "delete the entry " + clrtxt.bold + clrtxt.darkcyan + \
                      "'coordinates'" + clrtxt.end + " from the Thetis json " + \
                      "input file, or update the entry " + clrtxt.bold + \
                      clrtxt.darkcyan + "'number of CRS'" + clrtxt.end + " under" + \
                      " the entry " + clrtxt.bold + clrtxt.darkcyan + \
                      "'coordinates'" + clrtxt.end + " to the correct number of " + \
                      "available Coordinate Reference Systems"
            raise SystemExit(message)

        # If the number of available CRS is less than the number of coordinate
        # conversion reguested
        if no_crs <= no_conv:
            # Terminate Thetis
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "The number of available CRS is less than or " + \
                      "equal to the requested coordinate conversions. " + \
                      "Specifically, available CRS: " + clrtxt.bold + \
                      clrtxt.darkcyan + str(
                no_crs) + clrtxt.end + " , Number of" + \
                      " coordinate conversion requested: " + clrtxt.bold + \
                      clrtxt.darkcyan + str(no_conv) + clrtxt.end + ". The " + \
                      " number of available CRS should be one more than the " + \
                      "number of coordinate conversion requests"
            raise SystemExit(message)

        # If the number of available CRS is more than the number of coordinate
        # conversions requested
        elif no_crs - 1 > no_conv:
            # Terminate Thetis
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "The number of available CRS is more than the" + \
                      " requested coordinate conversions. Specifically, " + \
                      "available CRS: " + clrtxt.bold + clrtxt.darkcyan + \
                      str(no_crs) + clrtxt.end + " , Number of coordinate " + \
                      "conversion requested: " + clrtxt.bold + clrtxt.darkcyan + \
                      str(no_conv) + clrtxt.end + " The number of available CRS " + \
                      "should be one more than the number of coordinate " + \
                      "conversion requests"
            raise SystemExit(message)

        # Add mesh coordinates to the coordinates dictionary, which will cont-
        # ain all the available CRS for this simulation
        self.coord_dict = {
            self.config["mesh parameters"]["epsg"]: {
                "x-coords": self.xy[:, 0],
                "y-coords": self.xy[:, 1],
                "x-coord name": self.config["mesh parameters"]["x-coord name"],
                "y-coord name": self.config["mesh parameters"]["y-coord name"]
            }
        }

        print_output("- Convert mesh coordinates to other CRS")

        # Loop through the number of conversions to be performed
        for index in range(no_conv):
            # Check if the input CRS exists, i.e. we have the coordinates
            if dict[conv[index]]['input epsg'] in self.coord_dict:
                # Perform the coordinate conversion
                CoordsTools.process_coords_conversion(self, dict[conv[index]])
            # If the input CRS isn't available, i.e. it isn't in the
            # self.coord_dict
            else:
                # Terminate Thetis
                keys = ','.join(str(k) for k in list(self.coord_dict.keys()))
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                          clrtxt.end + "The input epsg, " + clrtxt.bold + \
                          clrtxt.darkcyan + str(
                    dict[conv[index]]['input epsg']) + \
                          clrtxt.end + ", to be used to convert coordinates to " + \
                          "the requested format does not exist in the already " + \
                          "available CRS: " + clrtxt.bold + clrtxt.darkcyan + \
                          str(keys) + clrtxt.end
                raise SystemExit(message)

        return


    def convert_coords_to_another_epsg(self, dict, input_epsg, output_epsg):
        """
        Convert coordinates from one CRS with a known EPSG to another CRS with
        a known EPSG utilising the pyproj.Transformer utility.
        Inputs:
        self        :
        dict        : A dictionary containing the necessary keys-values for co-
                      ordinate conversion from one EPS to another
        input_epsg  : An integer with the epsg code
        output_epsg : An integer with the epsg code

        Outputs:
        None

        Returns:
        self.coord_dict : Added another CRS in the self.coord_dict, which con-
                          tains all the available CRS and their coordinates.
        """
        from pyproj import CRS, Transformer

        # Specify the EPSG code for the coordinates we know (Input)
        inProj = CRS.from_epsg(input_epsg)
        # Specify the desired epsg code for the coordinates we want (Output)
        outProj = CRS.from_epsg(output_epsg)

        # Create a transformer object to convert the input coordinates to their de-
        # sired coordinate system
        transformer = Transformer.from_crs(
            crs_from=inProj, crs_to=outProj, always_xy=True
        )

        # Load known coordinates
        x_in_coords = self.coord_dict[input_epsg]["x-coords"]
        y_in_coords = self.coord_dict[input_epsg]["y-coords"]

        # Assign arrays for the unknown coordinates
        x_out_coords = np.zeros(np.shape(x_in_coords))
        y_out_coords = np.zeros(np.shape(y_in_coords))

        # Loop through the known coordinates:
        for index in range(np.shape(x_in_coords)[0]):
            # Get the truple with the unknown coordinates
            help = (transformer.transform(
                x_in_coords[index], y_in_coords[index])
            )
            # Pass the coordinates in the appropriate array
            x_out_coords[index] = help[0]
            y_out_coords[index] = help[1]

        # Add the new coordinates in the coordinates dictionary
        self.coord_dict[output_epsg] = {
            "x-coords": x_out_coords,
            "y-coords": y_out_coords,
            "x-coord name": dict['output x-coord name'],
            "y-coord name": dict["output x-coord name"]
        }

        return

    def convert_local_crs_to_latlon(self, dict, input_epsg, output_epsg):
        """
        Convert the coordinates of a local Cartesian coordinate System to a
        LatLon coordinate system utilising a rotation matrix
        | m_deg_long * (long - long_0) |   | cos(theta)  -sin(theta) |   | x |
        |                              | = |                         | * |   |
        | m_deg_lat * (lat - lat_0)    |   | sin(theta)   cos(theta) |   | y |
        Inputs:
        self        :
        dict        : A dictionary containing the necessary keys-values for co-
                      ordinate conversion from a local Cartesian system to Lat-
                      Lon
        input_epsg  : A string containing the name of the 'Local CRS'
        output_epsg : An integer with the epsg code

        Outputs:
        None

        Returns:
        self.coord_dict : Added another CRS in the self.coord_dict, which con-
                          tains all the available CRS and their coordinates.
        """
        ## Check that the required keys are within the dictionary
        # Load origin longitude if the entry exists
        if "origin longitude [o]" in dict:
            long_0 = dict["origin longitude [o]"]
        # If the appropriate key doesn't exist, terminate Thetis
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "No origin longitude is provided for the " + \
                      "conversion from " + clrtxt.bold + clrtxt.darkcyan + \
                      str(input_epsg) + clrtxt.end + " to the CRS with EPSG " + \
                      clrtxt.bold + clrtxt.darkcyan + str(
                output_epsg) + clrtxt.end
            raise SystemExit(message)

        # Load origin latitude, if the appropriate key exists
        if "origin latitude [o]" in dict:
            lat_0 = dict["origin latitude [o]"]
        # If the appropriate key doesn't exist, terminate Thetis
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "No origin latitude is provided for the " + \
                      "conversion from " + clrtxt.bold + clrtxt.darkcyan + \
                      str(input_epsg) + clrtxt.end + " to the CRS with EPSG " + \
                      clrtxt.bold + clrtxt.darkcyan + str(
                output_epsg) + clrtxt.end
            raise SystemExit(message)

        # Load meter per longitude degree, if the appropriate key exists
        if "meter per longitude degree" in dict:
            m_deg_long = dict["meter per longitude degree"]
        # If the appropriate key doesn't exist, terminate Thetis
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "No meters per longitude degree is provided " + \
                      "for the conversion from " + clrtxt.bold + clrtxt.darkcyan + \
                      str(input_epsg) + clrtxt.end + " to the CRS with EPSG " + \
                      clrtxt.bold + clrtxt.darkcyan + str(
                output_epsg) + clrtxt.end
            raise SystemExit(message)

        # Load meter per latitude degree if the appropriate key exists
        if "meter per latitude degree" in dict:
            m_deg_lat = dict["meter per latitude degree"]
        # If the appropriate key doesn't exist, terminate Thetis
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "No meters per latitude degree is provided " + \
                      "for the conversion from " + clrtxt.bold + clrtxt.darkcyan + \
                      str(input_epsg) + clrtxt.end + " to the CRS with EPSG " + \
                      clrtxt.bold + clrtxt.darkcyan + str(
                output_epsg) + clrtxt.end
            raise SystemExit(message)

        # Load angle from Local CRS to LatLon if the appropriate key exists
        if "angle from input to output epsg [o]" in dict:
            theta = dict["angle from input to output epsg [o]"]
            # Convert angle to radians
            theta = theta * np.pi / 180
        # If the appropriate key doesn't exist, terminate Thetis
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "No angle between the two CRS is provided " + \
                      "for the conversion from " + clrtxt.bold + clrtxt.darkcyan + \
                      str(input_epsg) + clrtxt.end + " to the CRS with EPSG " + \
                      clrtxt.bold + clrtxt.darkcyan + str(
                output_epsg) + clrtxt.end
            raise SystemExit(message)

        # Load local Cartesian coordinates from coordinates dictionary
        # Load local Cartesian coordinates from coordinates dictionary
        x = self.coord_dict[input_epsg]["x-coords"]
        y = self.coord_dict[input_epsg]["y-coords"]

        # Calculate the latlon coordinates
        long = (x * np.cos(theta) - y * np.sin(theta)) / m_deg_long + long_0
        lat = (x * np.sin(theta) + y * np.cos(theta)) / m_deg_lat + lat_0

        # Add coordinate to coordinate dictionary, for them to be available
        self.coord_dict[output_epsg] = {
            "x-coords": long,
            "y-coords": lat,
            "x-coord name": dict['output x-coord name'],
            "y-coord name": dict["output x-coord name"]
        }

        return


    def get_mesh_coordinates(self):
        """
        Get the mesh coordinates in (x,y) format and (y,x) format. The latter
        is for netCDF interpolation
        Inputs:
        self :

        Outputs:
        None

        Returns
        self.xy : 2-D array containing the mesh coordinates in (x,y) order
        self.yx : 2-D array containing the mesh coordinates in (y,x) order
        """
        # Get coordinates
        self.xy = self.mesh2d.coordinates.dat.data.copy()
        # Convert the order of the coordinates i.e. (y,x) for interpolation
        # from netcdf
        self.yx = self.xy[:, [1, 0]]

        return

    def get_mesh_coords_in_wave_char_crs(self, dict):
        """
        Get the mesh coordinates in the EPSG utilised in the netcdf containing
        the wave characteristics. If the EPSG is not available in the coordina-
        te dictionary, terminate Thetos
        Inputs:
        self:
        dict : A dictionary containing the coordinates parameters of the eta
               netcdf

        Outputs:
        None

        Returns:
        self.wave_coords : 2-D array containing the mesh coordinates in the
                           EPSG of the eta netcdf
        """
        # If the wind forcing utilise the same epsg as the mesh
        if dict["epsg"] == self.config["mesh parameters"]["epsg"]:
            self.wave_coords = self.yx

        # If the wind netcdf utilises a different than the mesh
        else:
            # If the coordinates are available in the coordinates dictionary
            if int(dict["epsg"]) in list(self.coord_dict.keys()):
                self.wave_coords = np.column_stack(
                    (
                        self.coord_dict[int(dict["epsg"])]["y-coords"],
                        self.coord_dict[int(dict["epsg"])]["y-coords"]
                    )
                )

            # If the coordinates aren't available, terminate Thetis
            else:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The EPSG" + clrtxt.bold + clrtxt.darkcyan + \
                    " " + str(dict["epsg"]) + clrtxt.end + " is not " + \
                    "available in the coordinates dictionary. Please add " + \
                    "the necessary coordinate conversion under the entry " + \
                    "' " + clrtxt.bold + clrtxt.darkcyan + "coordinates " + \
                    "conversion" + clrtxt.end + "'."
                raise SystemExit(message)

        return


    def get_mesh_coords_in_TPXO_crs(self, boundary_ids=None, TPXO_forcing=False):
        """
        Get the mesh coordinates in the EPSG 4326, i.e. in LatLon coordinates,
        which is the TPXO epsg. If the mesh epsg is not 4326, then conversion
        is needed.
        Inputs:
        self:

        Outputs:
        None

        Returns:
        self.tpxo_coords : 2-D array containing the mesh coordinates in EPSG
                           4326
        """
        # Get the epsg of the mesh coordinates
        mesh_epsg = self.config["mesh parameters"]["epsg"]
        # Specify the TPXO epsg
        tpxo_epsg = int(4326)

        # If boundary ids is not None, get the indices corresponfing to the
        # boundary nodes
        if boundary_ids is not None:
            # Get the number of owned elements, i.e. number of nodes
            no_owned = self.P1_2d.cell_node_map().toset.sizes[1]
            # Get the indices of nodes (included the ghosted ones, correspond-
            # ing to the boundary IDs
            self.bc_ind = DirichletBC(self.P1_2d, 0.0, boundary_ids).nodes
            # Find the indices of the ghosted nodes
            del_ind = np.where( no_owned <= self.bc_ind)
            if len(del_ind) != 0:
                self.bc_ind = np.delete(self.bc_ind, del_ind )



        # If the mesh epsg is a string, i.e. a local CRS
        if type(mesh_epsg) == str:
            # Get the available EPSG/CRS
            coords_epsg = list(self.coord_dict.keys())

            # If no coordinates in EPSG 4326 are available through coordinates
            # conversion, terminate Thetis
            if tpxo_epsg not in coords_epsg:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The mesh epsg is " + clrtxt.bold + \
                    clrtxt.darkcyan + "'Local CRS'" + clrtxt.end + " and " + \
                    "no coordinate conversion to " + clrtxt.bold + \
                    clrtxt.darkcyan + " EPSG 4326 " + clrtxt.end + "has " + \
                    "been specified for " + clrtxt.bold + clrtxt.darkcyan + \
                    "TPXO forcing" + clrtxt.end + ". Please add the " + \
                    "necessary CRS conversion under the entry '" + \
                    clrtxt.bold + clrtxt.darkcyan + "coordinates " + \
                    "conversion" + clrtxt.end + "'."
                raise SystemExit(message)

            # Get TPXO coordinates
            self.tpxo_coords = np.column_stack( (
                self.coord_dict[tpxo_epsg]["x-coords"][self.bc_ind],
                self.coord_dict[tpxo_epsg]["y-coords"][self.bc_ind]
            ))

        # If the EPSG is a number
        elif type(mesh_epsg) == int or type(mesh_epsg) == float:
            # If the mesh EPSG is not 4326
            if int(mesh_epsg) != tpxo_epsg:
                coords_epsg = list(self.coord_dict.keys())

                # Check that we have coordinates with EPSG 4326 available in
                # the coordinates dictionary, else terminate Thetis
                if tpxo_epsg not in coords_epsg:
                    message = clrtxt.bold + clrtxt.red + "Terminating " + \
                        "Thetis: " + clrtxt.end + "No coordinates with " + \
                        clrtxt.bold + clrtxt.darkcyan + "EPSG 4326 " + \
                        clrtxt.end + "are available. Please add the " + \
                        "necessary CRS conversion under the entry '" + \
                        clrtxt.bold + clrtxt.darkcyan + "coordinates " + \
                        "conversion" + clrtxt.end + "'."
                    raise SystemExit(message)

                # Get TPXO coordinates
                self.tpxo_coords = np.column_stack((
                    self.coord_dict[tpxo_epsg]["x-coords"][self.bc_ind],
                    self.coord_dict[tpxo_epsg]["y-coords"][self.bc_ind]
                ))

            # If the mesh EPSG is 4326
            else:
                # Pass the mesh coordinates to the TPXO ones.
                self.tpxo_coords = self.xy[self.bc_ind, :]

        # If TPXO forcing, we need to change the order of the coordinates
        if TPXO_forcing:
            self.tpxo_coords = np.column_stack((
                self.tpxo_coords[:,1],
                self.tpxo_coords[:,0]
            ))

        return

    def get_mesh_coords_in_uv_crs(self, dict):
        """
        Get the mesh coordinates in the EPSG utilised in the netcdf containing
        the uv forcing. If the EPSG is not available in the coordinate dict-
        ionary, terminate Thetos
        Inputs:
        self:
        dict : A dictionary containing the coordinates parameters of the uv
               netcdf

        Outputs:
        None

        Returns:
        self.wind_coords : 2-D array containing the mesh coordinates in the
                           EPSG of the uv netcdf
        """
        # If boundary ids is not None, get the indices corresponfing to the
        # boundary nodes
        if "boundary ids" in dict:
            # Get the boundary ids
            boundary_ids = dict["boundary ids"]
            # Get the number of owned elements, i.e. number of nodes
            no_owned = self.P1_2d.cell_node_map().toset.sizes[1]
            # Get the indices of nodes (included the ghosted ones, correspond-
            # ing to the boundary IDs
            self.bc_ind = DirichletBC(self.P1_2d, 0.0, boundary_ids).nodes
            # Find the indices of the ghosted nodes
            del_ind = np.where(no_owned <= self.bc_ind)
            if len(del_ind) != 0:
                self.bc_ind = np.delete(self.bc_ind, del_ind)
        else:
            # Get the number of owned elements, i.e. number of nodes
            no_owned = self.P1_2d.cell_node_map().toset.sizes[1]
            # Get all the nodes
            self.bc_ind = np.arange(no_owned)


        # If the wind forcing utilise the same epsg as the mesh
        if dict["epsg"] == self.config["mesh parameters"]["epsg"]:
            self.uv_coords = self.yx[self.bc_ind, :]

        # If the wind netcdf utilises a different than the mesh
        else:
            # If the coordinates are available in the coordinates dictionary
            if int(dict["epsg"]) in list(self.coord_dict.keys()):
                mes="Check this out the velocity coordinates out"
                raise SystemExit(message)
                self.uv_coords = np.column_stack(
                    (
                        self.coord_dict[int(dict["epsg"])]["y-coords"][self.bc_ind],
                        self.coord_dict[int(dict["epsg"])]["x-coords"][self.bc_ind]
                    )
                )

            # If the coordinates aren't available, terminate Thetis
            else:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The EPSG" + clrtxt.bold + clrtxt.darkcyan + \
                    " " + str(dict["epsg"]) + clrtxt.end + " is not " + \
                    "available in the coordinates dictionary. Please add " + \
                    "the necessary coordinate conversion under the entry " + \
                    "' " + clrtxt.bold + clrtxt.darkcyan + "coordinates " + \
                    "conversion" + clrtxt.end + "'."
                raise SystemExit(message)


        return


    def get_mesh_coords_in_wind_crs(self, dict):
        """
        Get the mesh coordinates in the EPSG utilised in the netcdf containing
        the wind forcing. If the EPSG is not available in the coordinate dict-
        ionary, terminate Thetos
        Inputs:
        self:
        dict : A dictionary containing the coordinates parameters of the wind
               netcdf

        Outputs:
        None

        Returns:
        self.wind_coords : 2-D array containing the mesh coordinates in the
                           EPSG of the wind netcdf
        """
        # If the wind forcing utilise the same epsg as the mesh
        if dict["epsg"] == self.config["mesh parameters"]["epsg"]:
            self.wind_coords = self.yx

        # If the wind netcdf utilises a different than the mesh
        else:
            # If the coordinates are available in the coordinates dictionary
            if int(dict["epsg"]) in list(self.coord_dict.keys()):
                self.wind_coords = np.column_stack(
                    (
                        self.coord_dict[int(dict["epsg"])]["y-coords"],
                        self.coord_dict[int(dict["epsg"])]["x-coords"]
                    )
                )

            # If the coordinates aren't available, terminate Thetis
            else:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The EPSG" + clrtxt.bold + clrtxt.darkcyan + \
                    " " + str(dict["epsg"]) + clrtxt.end + " is not " + \
                    "available in the coordinates dictionary. Please add " + \
                    "the necessary coordinate conversion under the entry " + \
                    "' " + clrtxt.bold + clrtxt.darkcyan + "coordinates " + \
                    "conversion" + clrtxt.end + "'."
                raise SystemExit(message)


        return


    def get_mesh_coords_in_swan_crs(self, dict):
        """
        Get the mesh coordinates in the EPSG utilised in the netcdf containing
        the wind forcing. If the EPSG is not available in the coordinate dict-
        ionary, terminate Thetos
        Inputs:
        self:
        dict : A dictionary containing the coordinates parameters of the wind
               netcdf

        Outputs:
        None

        Returns:
        self.wind_coords : 2-D array containing the mesh coordinates in the
                           EPSG of the wind netcdf
        """
        # If the wind forcing utilise the same epsg as the mesh
        if dict["epsg"] == self.config["mesh parameters"]["epsg"]:
            self.swan_coords = self.xy

        # If the wind netcdf utilises a different than the mesh
        else:
            # If the coordinates are available in the coordinates dictionary
            if int(dict["epsg"]) in list(self.coord_dict.keys()):
                self.swan_coords = np.column_stack(
                    (
                        self.coord_dict[int(dict["epsg"])]["x-coords"],
                        self.coord_dict[int(dict["epsg"])]["y-coords"]
                    )
                )

            # If the coordinates aren't available, terminate Thetis
            else:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The EPSG" + clrtxt.bold + clrtxt.darkcyan + \
                    " " + str(dict["epsg"]) + clrtxt.end + " is not " + \
                    "available in the coordinates dictionary. Please add " + \
                    "the necessary coordinate conversion under the entry " + \
                    "' " + clrtxt.bold + clrtxt.darkcyan + "coordinates " + \
                    "conversion" + clrtxt.end + "'."
                raise SystemExit(message)

        return


    def get_swan_coords_in_mesh_crs(self, dict):
        """
        Get the mesh coordinates in the EPSG utilised in the netcdf containing
        the wind forcing. If the EPSG is not available in the coordinate dict-
        ionary, terminate Thetos
        Inputs:
        self:
        dict : A dictionary containing the coordinates parameters of the wind
               netcdf

        Outputs:
        None

        Returns:
        self.wind_coords : 2-D array containing the mesh coordinates in the
                           EPSG of the wind netcdf
        """
        from pyproj import CRS, Transformer

        mesh_epsg = self.config["mesh parameters"]["epsg"]
        swan_epsg = dict["epsg"]
        # If the wind forcing utilise the same epsg as the mesh
        if dict["epsg"] == mesh_epsg:
            self.swan_coords_in_mesh_crs = self.xy_swan

        # If the wind netcdf utilises a different than the mesh
        else:
            # Specify the EPSG code for the coordinates we know (Input)
            inProj = CRS.from_epsg(swan_epsg)
            # Specify the desired epsg code for the coordinates we want (Output)
            outProj = CRS.from_epsg(mesh_epsg)

            # Create a transformer object to convert the input coordinates to their de-
            # sired coordinate system
            transformer = Transformer.from_crs(
                crs_from = inProj, crs_to = outProj, always_xy = True
            )

            x_in_coords = self.xy_swan[:, 0]
            y_in_coords = self.xy_swan[:, 1]

            # Assign arrays for the unknown coordinates
            x_out_coords = np.zeros(np.shape(x_in_coords))
            y_out_coords = np.zeros(np.shape(y_in_coords))

            # Loop through the known coordinates:
            for index in range(np.shape(x_in_coords)[0]):
                # Get the truple with the unknown coordinates
                help = (transformer.transform(
                    x_in_coords[index], y_in_coords[index])
                )
                # Pass the coordinates in the appropriate array
                x_out_coords[index] = help[0]
                y_out_coords[index] = help[1]

            self.swan_coords_in_mesh_crs = np.column_stack(
                (x_out_coords, y_out_coords)
            )


            # # If the coordinates aren't available, terminate Thetis
            # else:
            #     message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
            #         clrtxt.end + "The EPSG" + clrtxt.bold + clrtxt.darkcyan + \
            #         " " + str(dict["epsg"]) + clrtxt.end + " is not " + \
            #         "available in the coordinates dictionary. Please add " + \
            #         "the necessary coordinate conversion under the entry " + \
            #         "' " + clrtxt.bold + clrtxt.darkcyan + "coordinates " + \
            #         "conversion" + clrtxt.end + "'."
            #     raise SystemExit(message)
        return


    def process_coords_conversion(self, dict):
        """
        Determine which coordinate conversion is requested, check if it is fea-
        sible. If it is, call the appropriate method and update the coordinate
        dictionary to make the new CRS available to the whole BmiThetis class.
        Inputs:
        self :
        dict : A dictionary containing the parameters for the coordinate conve-
               sion requested.

        Outputs:
        None

        Returns:
        self.coord_dict : Added the requested CRS in the self.coord_dict, which
                          contains all the available CRS and their coordinates.
        """
        from pyproj import CRS

        # Get the input epsg
        input_epsg = dict["input epsg"]
        # Get the output epsg
        output_epsg = dict["output epsg"]

        print_output(f"  * Convert coordinates from {input_epsg} to" + \
                     f" {output_epsg}.")

        # If we are dealing with a local Cartesian Coordinate System as the in-
        # put CRS
        if input_epsg == 'Local CRS':
            # If the output coordinate system is also 'Local CRS', terminate
            # Thetis
            if output_epsg == 'Local CRS':
                message = clrtxt.bold + clrtxt.red + "Terminating " + \
                          "Thetis :" + clrtxt.end + "Conversion from a local" + \
                          " Cartesian System is only available to a LatLon " + \
                          "Coordinate System and not to another local Cartesian" + \
                          " coordinate system. Please provide an output epsg " + \
                          "code corresponding to LatLon Coordinates, such as" + \
                          " 4326"
                raise SystemExit(message)

            # Convert epsg code to CRS
            outProj = CRS.from_epsg(int(output_epsg))
            # If the output CRS is LatLon, convert/transform the coordinates
            if outProj.is_geographic:
                CoordsTools.convert_local_crs_to_latlon(
                    self, dict, input_epsg, output_epsg
                )
            # If the output CRS isn't LatLon, terminate Thetis
            else:
                message = clrtxt.bold + clrtxt.red + "Terminating " + \
                          "Thetis :" + clrtxt.end + "Conversion from a local" + \
                          " Cartesian System is only available to a LatLon " + \
                          "Coordinate System. Please provide an output epsg " + \
                          "code corresponding to LatLon Coordinates, such as" + \
                          " 4326"
                raise SystemExit(message)
            return

        # If the output coordinate system is a local Cartesian coordinate sys-
        # tem
        if output_epsg == "Local CRS":
            # If the input coordinate system is also 'Local CRS', terminate
            # Thetis
            if input_epsg == 'Local CRS':
                message = clrtxt.bold + clrtxt.red + "Terminating " + \
                          "Thetis :" + clrtxt.end + "Conversion to a local" + \
                          " Cartesian System is only available from a LatLon " + \
                          "Coordinate System & not from another local Cartesian" + \
                          " coordinate system. Please provide an output epsg " + \
                          "code corresponding to LatLon Coordinates, such as" + \
                          " 4326"
                raise SystemExit(message)

            # Convert epsg code to crs
            inProj = CRS.from_epsg(int(input_epsg))
            # If input CRS is LatLon, convert coordinates
            if inProj.is_geographic:
                # THIS NEEDS TO BE PROGRAMMED
                message = clrtxt.bold + clrtxt.red + "The conversion " + \
                          "from LatLon coordinates to a local Coordinate System " + \
                          "has not been programmed yet" + clrtxt.end
                raise SystemExit(message)
            # If the input CRS isn't LatLon, terminate Thetis
            else:
                message = clrtxt.bold + clrtxt.red + "Terminating " + \
                          "Thetis :" + clrtxt.end + "Conversion to a local" + \
                          " Cartesian System is only available from a LatLon " + \
                          "Coordinate System. Please provide an input epsg " + \
                          "code corresponding to LatLon Coordinates, such as" + \
                          " 4326"
                raise SystemExit(message)
            return

        # Convert coordinates from one epsg to another
        CoordsTools.convert_coords_to_another_epsg(
            self, dict, input_epsg, output_epsg
        )

        return


class BathymetryTools(BmiThetis):
    """
    Various utilities regarding bathymetry
    MEthods included in this class
    * bathymetry_field(self)
    * import_bathymetry_from_netCDF(self)
    * smoothen_bathymetry(self) : Recommended to avoid it
    """

    def bathymetry_field(self):
        """
        Determine through the options of the json file how to configure the ba-
        thymetry function and call the appropriate methods
        Inputs
        self :

        Outputs:
        None

        Returns:
        self.bathymetry2d : The bathymetry function for the simulation
        """
        # If the bathymetry is to be imported from a netcdf
        if self.config["bathymetry parameters"]["import bathymetry"]:
            # Import bathymetry
            BathymetryTools.import_bathymetry_from_netCDF(self)
        # If the bathymetry is not to be imported from a netCDF
        else:
            # Terminate Thetis, as nothing has been programmed
            message = clrtxt.bold + clrtxt.red + " Bathymetry not imported" + \
                      " from a netcdf has not been coded yet."
            raise SystemExit(message)

        return


    def import_bathymetry_from_netCDF(self):
        """
        Deefine a dunction for the bathymetry and then call the appropriate me-
        thod to interpolate the bathymetry to the mesh coordinates. Perform so-
        me checks to ensure the smooth and correct application of the bathyme-
        try interpolation, like making sure that both the mesh and the netCDF
        coordinates utilise the same CRS. If the user wants a csv with the int-
        erpolated bathymetry, create a csv.
        Inputs:
        self :

        Outputs:
        If the output options -> csv -> create is true, create a csv with the
        mesh coordinates and the interpolated bathymtery in the mesh coordina-
        tes

        Returns:
        self.bathymetry2d : The bathymetry function for the simulation
        """
        from pathlib import Path

        print_output("- Generate bathymetry field from NetCDF file")

        # Load dictionary regardin the bathymetry parameters
        dict = self.config["bathymetry parameters"]

        ## Bathymetry netCDF ##
        # Get bathymetric netcdf filepath
        bath_filepath = dict["nc file"]
        # If the specified netCDF doesn't exist, terminate Thetis
        if not Path(bath_filepath).exists():
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "The specified netCDF for the bathymetry " + \
                      clrtxt.bold + clrtxt.darkcyan + bath_filepath + clrtxt.end + \
                      " does not exist"
            raise SystemExit(message)

        ## Thetis bathymetry ##
        # Define function
        self.bathymetry2d = Function(self.P1_2d, name="bathymetry")
        # Function as editable field
        data = self.bathymetry2d.dat.data

        # Get the EPSG / CRS of the coordinate system employed in the netCDF
        nc_epsg = dict["epsg"]
        # If the EPSG isn't the string 'local CRS', convert the number to integer
        if nc_epsg != "Local CRS":
            nc_epsg = int(nc_epsg)

        # Check if the netCDF epsg is the same as the CRS of the mesh coordina-
        # tes
        if nc_epsg == self.config["mesh parameters"]["epsg"]:
            # The coordinates to be interpolated will be the mesh ones
            nc_coords = self.yx
        # If the CRS of the bathymetry netcdf is different than the mesh'es
        else:
            # Check if the coordinates of the netcdf are available
            if int(nc_epsg) in self.coord_dict:
                # Get the corresponding x- and y-coordinates
                nc_x_coords = self.coord_dict[nc_epsg]["x-coords"]
                nc_y_coords = self.coord_dict[nc_epsg]["y-coords"]

            # If the coordinates are not available in the coordinate dictionary
            else:
                # Terminate Thetis
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                          clrtxt.end + "The CRS of the bathymetry netCDF, " + \
                          clrtxt.bold + clrtxt.darkcyan + str(nc_epsg) + \
                          clrtxt.end + ", is not available. Please add the " + \
                          "necessary coordinates conversion under the key " + \
                          "'coordinates'"
                raise SystemExit(message)

            # Concatenate the coordinates with the y-coordinates in the 0th col
            nc_coords = np.column_stack((nc_y_coords, nc_x_coords))

        # If the grid of the netCDF is rotated
        if dict["rotated grid"]:
            # Interpolate the bathymetry from the netCDF to the coordinates
            data = IOTools.import_irregular_nc(
                self,
                bath_filepath,
                dict["x name"],
                dict["y name"],
                dict["z name"],
                dict["interpolation method"],
                nc_coords
            )
        # If the grid is regular structured
        else:
            # Create interpolator from netCDF
            intp = IOTools.import_regular_nc(
                self,
                bath_filepath,
                dict["x name"],
                dict["y name"],
                dict["z name"],
                dict["interpolation method"]
            )
            # Interpolate bathymetry
            data = intp(nc_coords)

        # Update function's values
        self.bathymetry2d.dat.data[:] = data

        # If the user requests bathymetry smoothing
        if "bathymetry smoothing" in dict:
            BathymetryTools.smoothen_bathymetry(self,
                                                dict["no of smoothening"])

        # Check if the user wants to create any output for the bathymetry field
        if "output options" in dict:
            IOTools.process_output_requests(
                self, dict["output options"], "bathymetry", self.bathymetry2d
            )

        return

    def smoothen_bathymetry(self, no_smooth):
        """
        Smoothen the bathymetry to minimise drastic changes in a short distance
        which are liable to create instabilities.
        However, given the coupling, any smoothing should happen during the
        creation of the bathymetry netCDF (see that you implement this in your
        scripts)
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.bathymetry2d : A smoothed bathymetry
        """
        for index in range(no_smooth):
            print_output('  * Apply bathymetry smoothing')

            # Define weakforms for solution of PDE
            v = TestFunction(self.P1_2d)

            # No idea what is going on
            massb = assemble(v * self.bathymetry2d * dx)
            massl = assemble(v * dx)

            # Compute the component-wise multiplication
            with massl.dat.vec as ml, massb.dat.vec as mb, \
                    self.bathymetry2d.dat.vec as sb:
                ml.reciprocal()
                sb.pointwiseMult(ml, mb)

        return



class IOTools(BmiThetis):
    """
    Various methods regarding Input/Output utilities aiding in Thetis simula-
    tion
    Methods included in this class:
    * export_csv(self, data, data_name, units, filepath)
    * export_h5(self, func, filepath)
    * export_pvd(self, func, filepath)
    * get_nc_vector_format_and_var_names(self, dict)
    * import_h5(self, input_dir, filename,
    * import_irregular_nc(
        self, filename, x_name, y_name, z_name, intp_method, xy
      )
    * import_regular_nc(self, filename, x_name, y_name, z_name, intp_method)
    * import_transient_irregular_vector_nc(
        self, filename, x_name, y_name, t_name, z1_name, z2_name, t_type
      )
    * import_transient_regular_vector_nc(
        self, filename, x_name, y_name, t_name, z1_name, z2_name, t_type
      )
    * load_json_file(self, filename)
    * process_output_requests(self, dict, field_name, data)
    * read_data_from_csv(self, filename, col, type)
    """

    def export_csv(self, data, data_name, units, filepath):
        """
        Create a csv containing the coordinates of the mesh and the values in
        these points of a variable/field.
        Inputs:
        self      :
        data      : 1-D array containing the data to be saved in the csv
        data_name : String cntaining the name of the data to be saved
        units     : String containing the units of the data to be saved
        filepath  : The filepath where the csv will be saved

        Outputs:
        A csv with the following format:
        column 0 : x-coordinates of the mesh following the mesh CRS
        column 1 : y-coordinates of the mesh following the mesh CRS
        column 2 : The values of the field in the aforementioned coordinates

        Returns:
        None
        """
        # Create dictionary for the csv
        dict = {
            self.config["mesh parameters"]["x-coord name"]: self.xy[:, 0],
            self.config["mesh parameters"]["y-coord name"]: self.xy[:, 1],
            data_name + " " + units: data
        }
        # Create DataFrame
        df = pd.DataFrame(dict)
        # Save to csv
        df.to_csv(filepath, index=None)

        return


    def export_h5(self, func, f_name, filepath):
        f"""
        Export a field in h5 format (probably to continue the simulation)
        Inputs:
        Inputs:
        self     :
        func     : The function to be exported
        f_name   : The name of the function to be exported
        filepath : The path where the hdf5 file is to be saved

        Outputs:
        A hdf5 of the provided function

        Returns:
        None
        """
        with CheckpointFile(filepath, "w") as afile:
            # Save the mesh
            afile.save_mesh(self.mesh2d)
            afile.save_function(func, name=f_name)
            afile.close()

        # # OPen the checkpoint file for writing, create a new one if it doesn't
        # # exist, or erase everything on it if it does exist
        # chk = CheckpointFile(filepath, 'w')# DumbCheckpoint(filepath, mode=FILE_CREATE)
        # # Save the mesh
        # chk.save_mesh(self.mesh2d, )
        # # STre the function in the file
        # chk.save_function(func, name=f_name)#chk.store(func, name=f_name)

        return


    def export_pvd(self, func, filepath):
        """
        Export a field in pvd format to view in paraview
        Inputs:
        self     :
        func     : The function to be exported
        filepath : The path where the paraview file is to be saved

        Outputs:
        A pvd (paraview file) of the provided function

        Returns:
        None
        """
        # Create pvd
        File(filepath).write(func)

        return


    def get_nc_vector_format_and_var_names(self, dict):
        """
        Determine whether the vector in the netCDF is described by its compon-
        enets or its magnitude and direction given in Cartesian convention
        Inputs:
        self :
        dict : A dictionary containing the vector information

        Outputs:
        None

        Returns:
        z1_label : a string containing the name of the 1st component describing
                   the vector
        z2_label : a string containing the name of the 1st component describing
                   the vector
        """
        # Get the format of the vector
        v_form = dict["format"]

        # If the vector is given in x- and y-components
        if v_form == "components":
            z1_label = dict["x-component label"]
            z2_label = dict["y-component label"]

        # If the vector is given as magnitude and directrion
        elif v_form == "magnitude & direction":
            z1_label = dict["magnitude label"]
            z2_label = dict["direction label [o]"]

        # If another format is given, terminate Thetis
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "The vector formats in a netCDF can be: " + \
                clrtxt.bold + clrtxt.darkcyan + "'components'" + clrtxt.end + \
                " and " + clrtxt.bold + clrtxt.darkcyan + "'magnitude & " + \
                "direction'" + clrtxt.end
            raise SystemExit(message)

        return z1_label, z2_label


    def import_h5(self, input_dir, filename, f_name, func):
        """

        """
        import pandas as pd
        with CheckpointFile(input_dir + filename, 'r') as afile:
            mesh = afile.load_mesh("firedrake_default")
            f = afile.load_function(mesh, name=f_name)

        # Pass the values to the actual functions
        func.dat.data[:] = f.dat.data[:]

        # if f_name != "elevation" or f_name != "velocity":
        #
        #     # Create dataframe
        #     dict ={
        #         "x" : mesh.coordinates.dat.data[:,0],
        #         "y" : mesh.coordinates.dat.data[:,1],
        #         "parameter": f.dat.data[:]
        #     }
        #
        #     df = pd.DataFrame(dict)
        #     df.to_csv(f"./processed/{f_name}_formhdf5.csv", index=False)
        #
        #     dict = {
        #         "x": self.mesh2d.coordinates.dat.data[:, 0],
        #         "y": self.mesh2d.coordinates.dat.data[:, 1],
        #         "parameter": func.dat.data[:]
        #     }
        #     df = pd.DataFrame(dict)
        #     df.to_csv(f"./processed/{f_name}_imported.csv", index=False)
        #     exit()


        return


    def import_irregular_nc(
            self, filename, x_name, y_name, z_name, intp_method, xy
    ):
        """
        Interpolate the z field from x- and y-coordinates to the provided xy
        coordinates using the specified interpolation method. The known coord-
        inates and the value of z-field in these coordinates are obtained from
        the specified netcdf
        Inputs:
        self        :
        filename    : A string containing the filepath of the imported netCDF
                      file
        x_name      : A string containing the name of the x-coordinates in the
                      netCDF file
        y_name      : A string containing the name of the y-coordinates in the
                      netCDF file
        z_name      : A string containing the name of the variable in the net-
                      CDF file
        intp_method : A string containing the interpolation method: 'linear',
                      'nearest' or 'cubic
        xy          : 2-D array containing the coordinates where the variable
                      is to be interpolated to

        Outputs:
        None

        Returns:
        data : 1-D array containing the interpolated values corresponding to xy
        """
        from netCDF4 import Dataset
        from scipy.interpolate import griddata

        print_output('  * Reading NetCDF file')

        # Check that we have all the records of the rotated grid
        if not self.config["bathymetry parameters"]["full coords record"]:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + " When the netCDF grid is rotated, all the " + \
                      "coordinates of the netCDF must be included, i.e. the " + \
                      clrtxt.bold + clrtxt.darkcyan + "'full coords record'" + \
                      clrtxt.end + " has to be " + clrtxt.bold + clrtxt.darkcyan + \
                      "true" + clrtxt.end
            raise SystemExit(message)

        # Open file
        nc = Dataset(filename)
        # Extract coordinates and data
        x_data = nc.variables[x_name][:].astype('float64').flatten()
        y_data = nc.variables[y_name][:].astype('float64').flatten()
        z_data = nc.variables[z_name][:].astype('float64').flatten()

        # Interpolate data
        data = griddata(
            (y_data, x_data),
            z_data,
            xy,
            method=intp_method
        )

        # Close file
        nc.close()

        return data


    def import_regular_nc(self, filename, x_name, y_name, z_name, intp_method):
        """
        Create an interpolation object/function from data read from a NetCDF
        file.

        Inputs:
        self        :
        filename    : A string containing the filepath of the imported netCDF
                      file
        x_name      : A string containing the name of the x-coordinates in the
                      netCDF file
        y_name      : A string containing the name of the y-coordinates in the
                      netCDF file
        z_name      : A string containing the name of the variable in the net-
                      CDF file
        intp_method : A string containing the interpolation method: 'linear',
                      'nearest' or 'cubic

        Outputs:
        None

        Returns:
        intp : The interpolation function created from the netCDF data
        """
        from netCDF4 import Dataset
        from scipy.interpolate import RegularGridInterpolator

        print_output('  * Reading NetCDF file')
        # Open file
        nc = Dataset(filename)
        # Extract coordinates and data
        x_data = nc.variables[x_name][:].astype('float64')
        y_data = nc.variables[y_name][:].astype('float64')
        z_data = nc.variables[z_name][:].astype('float64')

        # If the netCDF contains all the coordinates of the structured grid
        if self.config["bathymetry parameters"]["full coords record"]:
            # Get only the values which the structured grid is comprised
            x_data = x_data[0, :]
            y_data = y_data[:, 0]


        # Create interpolator
        intp = RegularGridInterpolator(
            (y_data, x_data),
            z_data,
            method=intp_method,
            fill_value=0.0,
            bounds_error=False
        )
        # Close file
        nc.close()

        return intp


    def import_transient_irregular_vector_nc(
            self, filename, x_name, y_name, t_name, z1_name, z2_name, t_type
    ):
        """

        """
        from netCDF4 import Dataset

        print_output('  * Reading NetCDF file')

        # Check that we have all the records of the rotated grid
        if not self.config["bathymetry parameters"]["full coords record"]:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + " When the netCDF grid is rotated, all the " + \
                "coordinates of the netCDF must be included, i.e. the " + \
                clrtxt.bold + clrtxt.darkcyan + "'full coords record'" + \
                clrtxt.end + " has to be " + clrtxt.bold + clrtxt.darkcyan + \
                "true" + clrtxt.end
            raise SystemExit(message)

        # Open file
        nc = Dataset(filename)
        # Extract coordinates
        x_data = nc.variables[x_name][:].astype('float64').flatten()
        y_data = nc.variables[y_name][:].astype('float64').flatten()
        # Extract time
        t_data = nc.variables[t_name][:].astype(t_type).flatten()
        # Extract vector information
        z1_data = nc.variables[z1_name][:].astype('float64').flatten()
        z2_data = nc.variables[z2_name][:].astype('float64').flatten()

        return x_data, y_data, t_data, z1_data, z2_data


    def import_transient_regular_vector_nc(
            self, filename, x_name, y_name, t_name, z1_name, z2_name, t_type,
            full_record
    ):
        """

        """
        from netCDF4 import Dataset

        print_output('  * Reading NetCDF file')

        # Check that we have all the records of the rotated grid
        # if not self.config["bathymetry parameters"]["full coords record"]:
        #     message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
        #         clrtxt.end + " When the netCDF grid is rotated, all the " + \
        #         "coordinates of the netCDF must be included, i.e. the " + \
        #         clrtxt.bold + clrtxt.darkcyan + "'full coords record'" + \
        #         clrtxt.end + " has to be " + clrtxt.bold + clrtxt.darkcyan + \
        #         "true" + clrtxt.end
        #     raise SystemExit(message)

        # Open file
        nc = Dataset(filename)
        # Extract coordinates
        x_data = nc.variables[x_name][:].astype('float64')
        y_data = nc.variables[y_name][:].astype('float64')
        # Extract time
        t_data = nc.variables[t_name][:].astype(t_type)
        # Extract vector information
        z1_data = nc.variables[z1_name][:].astype('float64')
        z2_data = nc.variables[z2_name][:].astype('float64')

        # If the netCDF contains all the coordinates of the structured grid
        if full_record:
            # Get only the values which the structured grid is comprised
            x_data = x_data[0, :]
            y_data = y_data[:, 0]

        return x_data, y_data, t_data, z1_data, z2_data

    def import_transient_regular_scalar_nc(
            self, filename, x_name, y_name, t_name, z_name, t_type,
            full_record
    ):
        """

        """
        from netCDF4 import Dataset

        print_output('  * Reading NetCDF file')

        # # Check that we have all the records of the rotated grid
        # if not self.config["bathymetry parameters"]["full coords record"]:
        #     message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
        #         clrtxt.end + " When the netCDF grid is rotated, all the " + \
        #         "coordinates of the netCDF must be included, i.e. the " + \
        #         clrtxt.bold + clrtxt.darkcyan + "'full coords record'" + \
        #         clrtxt.end + " has to be " + clrtxt.bold + clrtxt.darkcyan + \
        #         "true" + clrtxt.end
        #     raise SystemExit(message)

        # Open file
        nc = Dataset(filename)
        # Extract coordinates
        x_data = nc.variables[x_name][:].astype('float64')
        y_data = nc.variables[y_name][:].astype('float64')
        # Extract time
        t_data = nc.variables[t_name][:].astype(t_type)
        # Extract vector information
        z_data = nc.variables[z_name][:].astype('float64')

        # If the netCDF contains all the coordinates of the structured grid
        if full_record:
            # Get only the values which the structured grid is comprised
            x_data = x_data[0, :]
            y_data = y_data[:, 0]

        return x_data, y_data, t_data, z_data


    def load_json_file(self, filename):
        """
        Check if a json configuration file was specified as input and load it
        if it exists
        Inputs:
        self     :
        filename : The filepath of the Thetis json configuration file

        Outputs:
        None

        Returns:
        self.config : A dictionary containing the Thetis configuration parame-
                      ters
        """
        import json
        from pathlib import Path

        # If no configuration file is specified
        if filename is None:
            # Exit Thetis
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "No json input file for Thetis was specified." + \
                      "Please specify the path for the json input file inside " + \
                      "the parenthesis of " + clrtxt.bold + clrtxt.darkcyan + \
                      "initialize" + clrtxt.end
            raise SystemExit(message)
        # If a configuration file is specified
        else:
            # Check that the file exists
            if not Path(filename).exists():
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                          clrtxt.end + "The specified json file '" + clrtxt.bold + \
                          clrtxt.darkcyan + filename + clrtxt.end + "' does not " + \
                          "exist."
                raise SystemExit(message)

            # Load the json file
            with open(filename) as json_file:
                self.config = json.load(json_file)

        return


    def process_output_requests(self, dict, field_name, func):
        """
        Process the format of the output request and call the appropriate me-
        thod.
        Inputs:
        self       :
        dict       : A dictionary containing the output parameters for all the
                     output requests regarding this field
        field_name : A string containing the name of the field to be saved
        func       : Thetis function to be exported in various formats

        Outputs:
        None

        Return:
        None
        """
        from pathlib import Path

        # Get the format of the outputs
        keys = list(dict.keys())

        # Loop through the output requests
        for key in keys:
            print_output(f"  * Creating {key} for {field_name}")

            # Get the output firectory
            output_dir = dict[key]["output directory"]
            # Create directory if it doesn't already exist
            Path(output_dir).mkdir(exist_ok=True, parents=True)
            # Get the name of the output file
            filename = dict[key]["filename"]

            # If a csv is to be exported
            if key == 'csv':
                # Get the values of the field in the mesh coordinates
                data = func.dat.data[:]
                # Get the units of the field
                units = dict[key]["units"]

                IOTools.export_csv(
                    self, data, field_name, units, output_dir + filename
                )
            # If a paraview file is to be exported
            elif key == 'pvd':
                IOTools.export_pvd(
                    self, func, output_dir + filename
                )
            # For all the other export options
            else:
                # Terminated Thetis
                message = clrtxt.bold + clrtxt + "Terminating Thetis: " + \
                    "No output option " + clrtxt.bold + clrtxt.darkcyan + \
                    "has been programmed yet." + clrtxt.end
                raise SystemExit

        return

    def read_data_from_csv(self, filename, cols, data_type, rows=0):
        """

        """
        from pathlib import Path

        # Check that the filepath exists
        if not Path(filename).exists():
            mes = ct.br + "Terminating Thetis: " + ct.e + "The filepath " + \
                  ct.bdc + str(filename) + ct.e + "does not exist."
            raise SystemExit(mes)

        data = np.loadtxt(
            filename, delimiter=',', usecols=cols, skiprows=rows, dtype=data_type
        )

        return data




class SetupTools(BmiThetis):
    """
    Various methods to aid in setting up the Thetis simulation, like defining
    the horizontal eddy viscosity field, the manning coefficient field, etc
    The methods included in this class are:
    * beta_plane_approx(self)
    * configure_mesh(self)
    * coriolis_coords(self, dict)
    * coriolis_field(self)
    * coriolis_frequency(self)
    * create_irregular_transient_vector_interp(
        self, x_data, y_data, t_data, z1_data, z2_data
      )
    * create_regular_transient_vector_interp(
        self, x_data, y_data, t_data, z1_data, z2_data
      )
    * create_transient_vector_interpolator(
        self, x_data, y_data, t_data, z1_data, z2_data, grid_type
      )
    * define_function_spaces(self)
    * eikonal_equation(self, dict)
    * manning_field(self)
    * set_status(self)
    * swe_boundary_conditions(self)
    * user_defined_functions(self)
    * viscosity_field(self)
    * viscosity_sponge(self, dict)
    * viscosity_sponge_open_bnd(self, dict)
    """

    def beta_plane_approx(self):
        """
        Apply the beta-plane approximation, which accounts for the variation of
        the Coriolis parameter f along the latitude
        f = f_0 + \beta * y
        where f     : The coriolis parameter accounting for variation with lat-
                      itude [rad/s]
              f_0   : The Coriolis parameter at latitude \phi_0 [rad/s]
              \beta : The Rossby parameter which is
                      \beta = 2 * \Omega * cos(\phi) / R
                      where R is the mean radius of Earth [m]
             y      : The meridional distance from \phi_0, i.e the difference
                      in Northings [m]
        For more information:
        * https://en.wikipedia.org/wiki/Beta_plane
        * https://glossary.ametsoc.org/wiki/Beta-plane
        * https://glossary.ametsoc.org/wiki/Meridional
        * https://www.thefreedictionary.com/Meridional+distance#:~:text=
            the%20distance%20or%20departure%20from%20the%20meridian%3B%20the
            %20easting%20or%20westing.

        Inputs:
        self :

        Outputs:
        self.coriolis : Function containing the vraying with latitude coriolis
                        frequency
        """
        from pyproj import CRS, Transformer

        # Convert latitude to radians
        lat = np.radians(self.lat_0)

        # Calculate the Rossby parameter
        beta = 2 * self.earth_omega * cos(lat) / self.earth_r

        # Check if we have any CRS in Northing/Easting available
        keys = list(self.coord_dict.keys())
        # Boolean variable specifying if a UTM crs is available. True if it is
        # available; False if it isn't
        utm_epsg = False
        # Loop through the available CRS
        for key in keys:
            # If the key isn't string, i.e. isn't a EPSG code
            if not isinstance(key, str):
                # Convert EPSG to CRS
                Proj = CRS.from_epsg(int(key))
                # Check if the CRS isn't LatLong
                if Proj.is_projected:
                    # Retain the epsg code
                    out_epsg = int(key)
                    # Update value to True
                    utm_epsg = True
                    break

        # If there are no coordinates in UTM CRS available, i.e. Easting/North-
        # ing, terminate Thetis
        if not utm_epsg:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "There are no UTM coordinates available. " + \
                      "UTM coordinates are necessary for the beta-approximation " + \
                      "of the Coriolis parameter. Add the necessary conversion " + \
                      "under the entry " + clrtxt.bold + clrtxt.darkcyan + \
                      "'coordinates' " + clrtxt.end + "in the json file."
            raise SystemExit(message)

        # Convert the representative coordinates from LatLong to UTM
        transformer = Transformer.from_crs(
            crs_from=CRS.from_epsg(4326),
            crs_to=CRS.from_epsg(out_epsg),
            always_xy=True
        )
        east_0, north_0 = transformer.transform(self.long_0, self.lat_0)

        # Get the Northing mesh coordinates
        northings = self.coord_dict[out_epsg]["y-coords"]

        # Interpolate the beta-plane approximation of the Coriolis frequency in
        # the mesh
        data = self.cor_freq + beta * (northings - north_0)
        # Pass these values into the function
        self.coriolis.dat.data[:] = data
        

        return


    def configure_mesh(self):
        """
        Load the mesh from the provided msh file
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.mesh : Function(?) containing the mesh of the simulation
        """
        from pathlib import Path

        # Get mesh filepath from json file
        mesh_filepath = self.config["mesh parameters"]["mesh file"]
        # If the mesh file does not exist, terminate Thetis
        if not Path(mesh_filepath).exists():
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "The specified msh file " + clrtxt.bold + \
                      clrtxt.darkcyan + mesh_filepath + clrtxt.end + " does " + \
                      "not exist"
            raise SystemExit(messgae)
        # Load mesh
        self.mesh2d = Mesh(mesh_filepath)

        return


    def coriolis_coords(self, dict):
        """
        Get the latitude and longitude which represent the domain. The coordi-
        nates can be user-defined or be calculated as the mean of the maximum
        and minimum of the mesh coordinates
        Inputs:
        self :
        dict : A dictionary containing the coriolis parameters

        Outputs:
        None

        Returns:
        self.lat_0  : The 'representative' latitude of the domain
        self.long_0 : The 'representative' longitude of the domain

        """
        # Check if the latitude is user defined or we will calculate it from
        # the mesh coordinates
        if dict["latitude [o]"] == 'unknown':
            # Check that the WGS 84 (EPSG:4326) is available, else terminate
            # Thetis
            if 4326 not in self.coord_dict:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The EPSG" + clrtxt.bold + clrtxt.darkcyan + \
                    " 4326 " + clrtxt.end + " is not available in the " + \
                    "coordinates dictionary. Add the necessary coordinate " + \
                    "conversion under the entry " + clrtxt.bold + \
                    clrtxt.darkcyan + "'coordinates' " + clrtxt.end
                raise SystemExit(message)

            # Gather all the coordinates
            lat = self.comm.allgather(self.coord_dict[4326]["y-coords"])
            lat_gathered = np.concatenate(lat)

            # Calculate the centroid
            self.lat_0 = np.mean(lat_gathered)

        else:
            self.lat_0 = dict["latitude [o]"]


        self.long_0 = 0





        # # Check if the longitude is user defined
        # if dict["longitude [o]"] == "unknown":
        #     # Check that the WGS 84 (EPSG:4326) is available, else terminate
        #     # Thetis
        #     if 4326 not in self.coord_dict:
        #         message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
        #             clrtxt.end + "The EPSG" + clrtxt.bold + clrtxt.darkcyan + \
        #             " 4326 " + clrtxt.end + " is not available in the " + \
        #             "coordinates dictionary. Add the necessary coordinate " + \
        #             "conversion under the entry " + clrtxt.bold + \
        #             clrtxt.darkcyan + "'coordinates' " + clrtxt.end
        #         raise SystemExit(message)
        #
        #     # Calculate the mean longitude of the domain according to the mesh
        #     # coordinates
        #     self.long_0 = np.mean(
        #         [np.amin(self.coord_dict[4326]["x-coords"]),
        #         np.amax(self.coord_dict[4326]["x-coords"])]
        #     )
        # else:
        #     self.long_0 = dict["longitude [o]"]

        return


    def coriolis_field(self):
        """
        Create a field containing the Coriolis frequency varying with latitude.
        To do so get or calculate the representative coordinates of the domain
        in LatLon, calculate the coriolis frequency at said coordinates and ap-
        ply the beta-plane approximation to get the varying Coriolis frequency.
        If the user has requested to export the Coriolis field, call the appr-
        opiate method
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.coriolis : Function containing the vraying with latitude coriolis
                        frequency
        """
        print_output("  * Generating Coriolis forcing field")

        # Load dictionary with coriolis parameters
        dict = self.config["coriolis parameters"]

        # Create coriolis forcing function
        self.coriolis = Function(self.P1_2d, name="coriolis")

        # The Earth's rotation rate [rad/s]
        self.earth_omega = 7.2921 * 10 ** (-5)
        # The Earth's mean radius [m]
        self.earth_r = 6371e3

        # Get the "mean"/given/representative coordinates of the domain
        SetupTools.coriolis_coords(self, dict)

        # Calculate coriolis frequency (a.k.a. Coriolis parameter, or Coriolis
        # frequency)
        SetupTools.coriolis_frequency(self)

        # Apply beta-plane approximation for the Coriolis
        SetupTools.beta_plane_approx(self)

        # Check if the user requested output for the Coriolis frequency
        if "output options" in dict:
            IOTools.process_output_requests(
                self,
                dict["output options"],
                "coriolis_frequency",
                self.coriolis
            )

        return


    def coriolis_frequency(self):
        """
        Calculate the Coriolis frequency f, also called the Coriolis parameter
        or Coriolis coefficient, which is equal to twice the rotation rate
        \Omega of the earth multiplied by the sine of the latitude \phi
        f = 2 * \Omega * sin(\phi)
        where f      : The Coriolis frequency
              \Omega : The rotation rate of the Eart equal to 7.2921 * 10^(-5)
                       rad/s
              \phi   : The latitude in degrees
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.cor_freq : A float number depicting the Coriolis frequency f
        """
        # Convert latitude to radians
        lat = np.radians(self.lat_0)

        # Calculate the Coriolis frequency
        self.cor_freq = 2 * self.earth_omega * sin(lat)

        return


    def create_regular_transient_vector_interp(
            self, x_data, y_data, t_data, z1_data, z2_data
    ):
        """
        Create a linear interpolation in 2-D space (x, y) and in time t for the
        components of a vector. The grid of the vector is structured regular.
        Inputs
        self :
        x_data  :
        y_data  :
        t_data  :
        z1_data :
        z2_data :

        Outputs:
        None

        Returns:
        z1_int : Interpolator for the x-component of the vector
        z2_int : Interpolator for the y-component of the vector
        """
        from scipy.interpolate import RegularGridInterpolator

        # Create the interpolator for x-component
        z1_int = RegularGridInterpolator(
            (t_data, y_data, x_data),
            z1_data,
            method = 'linear',
            bounds_error = False,
            fill_value = 0.0
        )

        # Create the interpolator for the y-component
        z2_int = RegularGridInterpolator(
            (t_data, y_data, x_data),
            z2_data,
            method='linear',
            bounds_error=False,
            fill_value=0.0
        )

        return z1_int, z2_int


    def create_transient_vector_interpolator(
            self, x_data, y_data, t_data, z1_data, z2_data, grid_type
    ):
        """
        Depending on the grid type of the vector, call the appropriate method
        to create the interpolator of the vector components x,y, across the 2-D
        space (x, y) and in time (t)
        Inputs
        self :
        x_data    :
        y_data    :
        t_data    :
        z1_data   :
        z2_data   :
        grid_type :

        Outputs:
        None

        Returns:
        z1_int : Interpolator for the x-component of the vector
        z2_int : Interpolator for the y-component of the vector

        """
        # If the grid of the vector isn't regular:
        if grid_type:
            interp = SetupTools.create_irregular_transient_vector_interp(
                self, x_data, y_data, t_data, z1_data, z2_data
            )

        # Else if the grid is regular
        else:
            z1_int, z2_int = SetupTools.create_regular_transient_vector_interp(
                self, x_data, y_data, t_data, z1_data, z2_data
            )

        return z1_int, z2_int

    def create_transient_scalar_interpolator(
            self, x_data, y_data, t_data, z_data, grid_type
    ):
        """
        Depending on the grid type of the vector, call the appropriate method
        to create the interpolator of the vector components x,y, across the 2-D
        space (x, y) and in time (t)
        Inputs
        self :
        x_data    :
        y_data    :
        t_data    :
        z1_data   :
        z2_data   :
        grid_type :

        Outputs:
        None

        Returns:
        z1_int : Interpolator for the x-component of the vector
        z2_int : Interpolator for the y-component of the vector

        """
        # If the grid of the vector isn't regular:
        if grid_type:
            print("THIS needs to be coded")
            # interp = SetupTools.create_irregular_transient_scalar_interp(
            #     self, x_data, y_data, t_data, z_data
            # )

        # Else if the grid is regular
        else:
            z_int = SetupTools.create_regular_transient_scalar_interp(
                self, x_data, y_data, t_data, z_data
            )

        return z_int

    def create_regular_transient_scalar_interp(
            self, x_data, y_data, t_data, z_data
    ):
        """
        Create a linear interpolation in 2-D space (x, y) and in time t for the
        components of a vector. The grid of the vector is structured regular.
        Inputs
        self :
        x_data  :
        y_data  :
        t_data  :
        z1_data :
        z2_data :

        Outputs:
        None

        Returns:
        z1_int : Interpolator for the x-component of the vector
        z2_int : Interpolator for the y-component of the vector
        """
        from scipy.interpolate import RegularGridInterpolator

        # Create the interpolator for x-component
        z_int = RegularGridInterpolator(
            (t_data, y_data, x_data),
            z_data,
            method='linear',
            bounds_error=False,
            fill_value=0.0
        )


        return z_int


    def define_function_spaces(self):
        """
        Define the functionspaces that may be used in the Thetis simulation.
        These are the Continuous Galerkin function-space for scalar, the Disco-
        ntinuous Galerking function-space for scalar, the Continous Galerkin
        function-space for vector; and the Discontinous Galerkin function-space
        for vector
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.P1DG_2d  : Discontinuous Galerkin function space for scalars
        self.P1_2d    : Continuous Galerkin function space for scalar
        self.P1v_2d   : Discontinuous Galerkin function space for vectors
        self.P1DGv_2d : Continuous Galerkin function space for vectors
        """
        # Define function spaces
        self.P1DG_2d = get_functionspace(self.mesh2d, "DG", 1)
        self.P1_2d = get_functionspace(self.mesh2d, "CG", 1)
        self.P1v_2d = VectorFunctionSpace(self.mesh2d, "CG", 1)
        self.P1DGv_2d = VectorFunctionSpace(self.mesh2d, "DG", 1)

        return


    def eikonal_equation(self, dict):
        """
        Solve the eikonal equation to determine the shortest path the wavefront
        has to travel in our domain.
        For more information on the eikonal equation:
        * https://en.wikipedia.org/wiki/Eikonal_equation
        * Slawinski, M.A., 2003. Seismic waves and rays in elastic media. Else-
        vier. Chapters 7-8
        Check if the user has requested to export the solution of the eikonal
        equation and call the appropriate method
        Inputs:
        self :
        dict : A dictionary containing the various parameters to be used in so-
               lving the Eikonal equation

        Outputs:
        None

        Returns:
        self.dist_bnd : A function containing the solution of the Eikonal equa-
                        tion depicting the shorted distance
        """
        print_output("  * Solving the Eikonal equation")

        # Check that we have 'ocean boundary' under PhysIDs
        if not "ocean boundary" in self.config["mesh parameters"]["PhysIDs"]:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "No PhysID with the name " + clrtxt.bold + \
                      clrtxt.darkcyan + 'ocean boundary' + clrtxt.end + " is " + \
                      "included under the entry " + clrtxt.bold + \
                      clrtxt.darkcyan + "'mesh parameters' " + clrtxt.end + \
                      "although the eikoanl equation is to be used for the " + \
                      "viscosity sponge"
            raise SystemExit(message)

        # Implement strong Dirichlet boundary condition in the ocean boundary
        bcs = [ DirichletBC(self.P1_2d, 0.0, dict["boundary ids"] ) ]

        # Create weak-forms for the solution of PDE
        self.dist_bnd_test = TestFunction(self.P1_2d)
        self.dist_bnd = Function(self.P1_2d)

        # Solve a Laplace equation for an initial guess of the solution of the
        # Eikonal equation
        F = dict["typical length"] ** 2 \
            * (inner(grad(self.dist_bnd), grad(self.dist_bnd_test))) * dx \
            - self.dist_bnd_test * dx
        # Solve the above equation
        solve(
            F == 0,
            self.dist_bnd,
            bcs,
            solver_parameters=dict["Laplace solver parameters"]
        )

        # Loop through eps values which determine the accuracy of the solution
        for eps in dict["eps"]:
            print_output(f"    - Solving the Eikonal equation with eps = " + \
                         f"{eps:.2f}")
            F = inner(
                sqrt(inner(grad(self.dist_bnd), grad(self.dist_bnd))),
                self.dist_bnd_test
            ) * dx - self.dist_bnd_test * dx + \
                eps * inner(grad(self.dist_bnd), grad(self.dist_bnd_test)) * dx
            solve(
                F == 0,
                self.dist_bnd,
                bcs,
                solver_parameters=dict["Eikonal solver parameters"]
            )

        # Check if the user requested output for the eikonal equation
        if "output options" in dict:
            if "eikonal equation" in dict["output options"]:
                IOTools.process_output_requests(
                    self,
                    dict["output options"]["eikonal equation"],
                    "eikonal_equation",
                    self.dist_bnd
                )

        return


    def manning_field(self):
        """
        Determine through the options of the json file if the manning coeffi-
        cient is constant across the domain or varies. If the former assign the
        value to the function, defined here, otherwise call the appropriate
        method. Check if the user has requested to export the manning coeffic-
        ient
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.n_mann : A function containing the manning coefficient field
        """
        # Load the dictionary with the manning parameters
        dict = self.config["manning parameters"]

        # Create manning function
        self.n_mann = Function(self.P1_2d, name="manning")

        # If the manning coefficient varies across the domain
        if dict["varying manning"]:
            # Terminate Thetis, as this has not been programmed yet
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: A " + \
                "varying manning coefficient field has not been programmed" + \
                "yet" + clrtxt.end
            raise SystemExit(message)

        # If the manning coefficient is constant across the domain
        else:
            print_output("- Generate constant manning field")
            # Assign the constant value to the function
            self.n_mann.assign(dict["value"])

        # Check if manning is to be exported:
        if "output options" in dict:
            IOTools.process_output_requests(
                self,
                dict["output options"],
                "manning_coeff",
                self.n_mann
            )

        return


    def set_status(self):
        """
        Deterine the coupling sttaus and the simulation kind

        Returns:
        self.simul_kind : A string containing the simulation kind
        """
        from bmi_tools import WCITools
        # Get the simulation kind
        self.simul_kind = self.config["simulation kind"]

        # Get the coupling status
        if "WCI" in self.config:
            WCITools.get_coupling_status(self)

        else:
            # Set the coupling status to "no_coupling
            self.coupl_stat = "no_coupling"


        # Set the ramp_coeff
        if self.simul_kind == "ramp" or self.simul_kind == "ramp_cont":
            # Get the duration of the ramping gradient
            self.ramp_grad = \
                self.config["ramp parameters"]["gradient duration [s]"]
            self.ramp_t_start = self.config["ramp parameters"]["start ramp time [s]"]

        # Print ramp coeff
        ThetisTools.calculate_ramp_coeff(self, 0)
        # if self.simul_kind == "ramp":
        #     self.ramp_coeff = 0.0
        # elif: self.simul_kind == "ramp_continue":
        #     self.ramp_coeff = 53.0
        # else:
        #     self.ramp_coeff = 1.0


        return



    def swe_boundary_conditions(self):
        """
        Create a dictionay for the Boundary conditions regarding the  Shallow-
        Water Equations based on the 'shallow-water BCs' entry in the json file.
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.swe_bnd : A dictionary contanining theboundary conditions regard-
                       ing the Shallow-water equations to be eployed in this
                       simulation
        """
        # Load dictionary containing the Shallow-water Boundary conditions
        dict = self.config["shallow-water BCs"]

        print_output("- Assigning boundary conditions")

        # Initialise dictionary
        self.swe_bnd = {}
        # Names of boundaries and boundary conditions
        bnd_names = list(dict.keys())
        bnd_cond = list(dict.values())

        # Correspondence between boundary names and PhysIDs
        phys_ids = self.config["mesh parameters"]["PhysIDs"]

        # Loop though the boundary names
        for index in range(np.shape(bnd_names)[0]):
            # If the boundary name doesn't exist in the Phys_id dictionary, te-
            # rminate Thetis
            if bnd_names[index] not in phys_ids:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The boundary name " + clrtxt.bold + \
                    clrtxt.darkcyan + f"'{bnd_names[index]}' " + clrtxt.end + \
                    "is not included in the entry " + clrtxt.bold + \
                    clrtxt.darkcyan + "'PhysIDs' " + clrtxt.end + "under " + \
                    "the entry " + clrtxt.bold + clrtxt.darkcyan + "'mesh" + \
                    " parameters' " + clrtxt.end
                raise SystemExit(message)

            # Initialise boundary condition dictionary for this particular bou-
            # ndary ID
            bc_dict = {}

            for no_bnd in range(np.shape(list(bnd_cond[index].items()))[0]):
                # Get boundary condition type and value
                bc_type = list(bnd_cond[index].items())[no_bnd][0]
                bc_value = list(bnd_cond[index].items())[no_bnd][1]

                # If the boundary condition value is a string, then we assign
                # a function (probably user-defined) in that boundary
                if isinstance(bc_value, str):
                    # Get the function and assign it in the boundary
                    bc_dict[bc_type] = getattr(self, bc_value)

                # If the boundary condition is a number
                else:
                    bc_dict[bc_type] = Constant(bc_value)

            # Update the Shallow-water BCs dictionary for the boundary ID
            self.swe_bnd[phys_ids[bnd_names[index]]] = bc_dict

        return


    def user_defined_function(self):
        """
        Define functions specified by the user under the entry 'function defin-
        ition' in the json file to be used in the simulation
        Inputs:
        self :

        Outputs:
        None

        Returns:
        The user-defined functions
        """
        # Load dictionary with the function definition parameters
        dict = self.config["function definition"]

        print_output("- Creating user-defined functions")

        # User-defined function names
        f_names = list(dict.keys())

        # Loop through the functions to be defined
        for name in f_names:
            print_output("  * Creating function " + name)

            # Get dictionary describing the function
            f_dict = dict[name]

            # Create the function
            setattr(
                self,
                name,
                Function(
                    getattr(  self,
                        f_dict["function_space"]
                    ),
                    name=f_dict["name"]
                )
            )

        return


    def viscosity_field(self):
        """
        Determine through the options of the json file if the viscosity is con-
        stant across the domain or if a viscosity sponge should be implemented
        and call the appropriate method or assign the constant value. The defi-
        nition of the horizontal eddy viscosity function takes place here.
        Also check if the user has requested to export the viscosity.
        Inputs:
        self

        Outputs:
        None

        Returns:
        self.h_visc : A function containing the horizontal eddy viscosity field
        """
        # Load the dict with the viscosity sponge parameters
        dict = self.config["viscosity parameters"]

        # Create viscosity function
        self.h_visc = Function(self.P1_2d, name="viscosity")

        # If a viscosity sponge should be implemented
        if dict["viscosity sponge"]:
            SetupTools.viscosity_sponge(self, dict)

        # If the eddy viscosity is constant across the domain
        else:
            print_output("- Generate constant viscosity field")

            if "value" not in dict:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The eddy viscosity is constant across " + \
                    "the domain but the key '" + clrtxt.bold + \
                    clrtxt.darkcyan + "value" + clrtxt.end + "' is not " + \
                    "included under the key '" + clrtxt.bold + \
                    clrtxt.darkcyan + "viscosity parameters" + clrtxt.end + \
                    "'."
                raise SystemExit(message)

            # Assign function's values
            self.h_visc.assign(dict["value"])

        # Check if viscosity is to be exported
        if "output options" in dict:
            if "viscosity" in dict["output options"]:
                IOTools.process_output_requests(
                    self,
                    dict["output options"]["viscosity"],
                    "viscosity",
                    self.h_visc
                )

        return


    def viscosity_sponge(self, dict):
        """
        Determine the 'format' of the viscosity sponge and call the appropriate
        method.
        Inputs:
        self :
        dict : A dictionary containing the parameters for the creation of a vi-
               scosity sponge

        Outputs:
        None

        Returns:
        None
        """
        print_output("- Creating a horizontal eddy viscosity sponge")

        # Check that the entry "eikonal equation" exists in the provided dict-
        # ionary
        if "use eikonal equation" not in dict:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "The key '" + clrtxt.bold + clrtxt.darkcyan + \
                      "use eikonal equation" + clrtxt.end + "' is not included " + \
                      "under the entry '" + clrtxt.bold + clrtxt.darkcyan + \
                      "viscosity parameters" + clrtxt.end + "', although the " + \
                      "option " + clrtxt.bold + clrtxt.darkcyan + "viscosity " + \
                      "sponge " + clrtxt.end + "is " + clrtxt.bold + \
                      clrtxt.darkcyan + "true" + clrtxt.end
            raise SystemExit(message)

        # ADD: Check that no entry 'value' exists

        # If the eikonal equation will be used for the creation of the viscosi-
        # ty sponge
        if dict["use eikonal equation"]:
            # Solve the Eikonal Equations
            SetupTools.eikonal_equation(self, dict)
            # Create a viscosity sponge alongside the open (ocean) boundaries
            # of the domain
            SetupTools.viscosity_sponge_open_bnd(self, dict)
        # If the eiknaol equation isn't used for the creation of a viscosity
        # sponge
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                "A viscosity sponge without using the eikonal equation has" + \
                " not been programmed yet."
            raise SystemExit(message)

        ## This is from Roelvink case where a user defined viscosity sponge was
        ## implemented
        # if self.config["viscosity parameters"]["where"] == "left":
        #     # Create function for interpolation
        #     dist_fn = Function(self.P1_2d, name="distance")
        #     dist_data = dist_fn.dat.data
        #     x_coords = self.xy[:, 0]
        #     y0 = self.config["viscosity parameters"]["max_value"]
        #     y1 = self.config["viscosity parameters"]["main_value"]
        #     dist = self.config["viscosity parameters"]["distance"]
        #     x0 = self.config["viscosity parameters"]["x0"]
        #     x1 = x0 + dist
        #     for index in range(np.shape(x_coords)[0]):
        #         if x_coords[index] > x1:
        #             dist_data[index] = 0
        #         else:
        #             dist_data[index] = 1 - (x_coords[index] - x0) / dist
        #     dist_fn.dat.data[:] = dist_data
        #
        #     print_output("- Generating vscosity sponge on the left side")
        #
        #     self.h_visc.interpolate(
        #         Max(
        #             y1, Constant(y0) * dist_fn
        #         )
        #     )
        #
        # IOTools.export_pvd(self, "viscosity", self.h_visc)

        return

    def viscosity_sponge_open_bnd(self, dict):
        """
        Create a viscosity sponge at the ocean boundaries utilising the solut-
        ion of the eikonal equation.
        Inputs:
        self :
        dict : A dictionary containing the parameters for the creation of a vi-
               scosity sponge

        Outputs:
        None

        Returns:
        self.h_visc : The function containing the eddy viscosity sponge
        """
        print_output("  * Creating viscosity sponge field")

        # Grab the maximum value of the viscosity sponge, the 'actual' value of
        # the viscosity and the distance where the values change
        main_value = dict["main value"]
        max_value = dict["max value"]
        dist = dict["gradient distance"]

        # Interpolate the viscosity values
        self.h_visc.interpolate(
            Max(
                main_value,
                max_value * (1 - self.dist_bnd / dist)
            )
        )

        return



class ThetisTools(BmiThetis):
    """
    Various methods for the configuration of the Flowsolver,
    This class contains the following methods:
    * calculate_ramp_coeff(self)
    * check_initial_conditions_names(self, names)
    * check_solver_times(self, dict)
    * create_flowsolver(self) under-constructions
    * get_export_params(self)
    * get_initial_conditions(self)
    * get_start_and_end_times(self, dict)
    """

    def __init__(self):
        """

        """
        BmiThetis._ic_names_options = [
            "elevation_ic",
            "u_velocity_ic",
            "v_velocity_ic",
            "uv_velocities_ic"
        ]
        BmiThetis._ic_wave_names = [
            "hwave_old",
            "dir_old",
            "qb_old",
            "wlen_old"
        ]

    def calculate_ramp_coeff(self, t):
        """

        """
        # Specify the ramp coefficient
        # If the simulation kind is ramp
        if (
                ((self.simul_kind == "ramp") and (t <= self.ramp_grad )) or
                ((self.simul_kind == "ramp_cont") and
                 (t + self.ramp_t_start <= self.ramp_grad ))
        ):
            self.ramp_coeff = \
                np.tanh((t + self.ramp_t_start) / (self.ramp_grad / 6.5))


        # If the simulation kind is 'new'
        elif (
                (self.simul_kind == "new") or
                (self.simul_kind == "continue") or
                ((self.simul_kind == "ramp") and (t > self.ramp_grad)) or
                ((self.simul_kind == "ramp_cont") and
                 (t + self.ramp_t_start > self.ramp_grad))
        ):
            self.ramp_coeff = 1.0



        return


    def check_initial_conditions_names(self, names):
        """
        Check that the keys under the entry "initial conditions" in the Thetis
        json configuration file are correct if the ICs are user-defined. The
        appropriate names are : elevation_ic, u_velocity_ic, v_velocity_ic. If
        one of them does not exist, terminate Thetis
        Inputs:
        self :
        names : A list containing the names of the initial condition functions
                as defined by the user

        Outputs:
        None

        Returns:
        None
        """
        # Check if we have coupling with Thetis
        if (
                ((self.coupl_stat=="2-way") or
                 (self.coupl_stat=="SWAN-to-Thetis")
                ) and
                ((self.simul_kind=="continue") or
                 (self.simul_kind=="ramp_cont"))
        ):
            # Remove the wave ic names
            for name in self._ic_wave_names:
                names.remove(name)

        # Loop through the names
        for name in names:
            # If the name is not one of the appropriate ones
            if name not in self._ic_names_options:
                mes = ct.br + "Terminating Thetis: " + ct.e + "The " + \
                      "function name in Thetis, i.e. the key containing " + \
                      "the information for the IC function, " + ct.bdc + \
                      str(name) + ct.e + " is not a viable option. Please " + \
                      "choose one of the following: " + ct.bdc + \
                      str(self._ic_names_options)[1:-1] + ct.e + "."
                raise SystemExit(mes)

        ## Check for WCI?
        return


    def check_solver_times(self, dict):
        """
        Perform some basic checks regarding the times of the simulation, i.e.
        if the timestep of the simulation is a multiple of the simulation dur-
        ation, if the export timestep is at least equal or bigger than the sim-
        ulation timestep and if the export timestep is a multiple of the simul-
        ation timestep.
        Inputs:
        self :
        dict : A dictionary containing the flowsolver parameters

        Outputs:
        None

        Returns:
        None
        """
        # Get the simulation timestep
        dt = dict["timestep"]
        exp_dt = dict["simulation_export_time"]

        # Check that the timestep is a factor of the simulation duration
        if self.t_duration % dt != 0:
            # Calculate the remainder
            rem = self.t_duration % dt
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "The simulation timestep " + clrtxt.bold + \
                clrtxt.darkcyan + str(dt) + clrtxt.end + " s is not a " + \
                "factor of the simulation duration " + clrtxt.bold + \
                clrtxt.darkcyan + str(self.t_duration) + clrtxt.end + " s." + \
                "The remainder is " + clrtxt.bold + clrtxt.darkcyan + \
                str(rem) + clrtxt.end + " s."
            raise SystemExit(message)

        # Check that the the export timemestep is or equal to the simulation
        # timestep
        if exp_dt < dt:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "The simulation timestep " + clrtxt.bold + \
                clrtxt.darkcyan + str(dt) + clrtxt.end + " s should not " + \
                "be lesser than the export timestep " + clrtxt.bold + \
                clrtxt.darkcyan + str(exp_dt) + clrtxt.end + " s."
            raise SystemExit(message)

        if exp_dt % dt != 0:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "The export timestep " + clrtxt.bold + \
                clrtxt.darkcyan + str(exp_dt) + clrtxt.end + " s should " + \
                "be a multiple of the simulation timestep " + clrtxt.bold + \
                clrtxt.darkcyan + str(dt) + clrtxt.end + " s."
            raise SystemExit(message)

        return


    def create_flowsolver(self):
        """
        Create the flowsolver object and assign all the necessary options based
        on the json file. This means, if the user has defined:
        * Define the flow solver
        * Calculate the start, end times and the duration of the simulation
        * Add viscosity field
        * Add wetting and drying field
        * Add manning coefficient field
        * Update the flowsolver options
        * Update the Shallow-Water Equation timestepper options
        * Assign the boundary conditions
        * Assign the initial conditions
        * Assign a detectors callback
        * Create a WCI object to account for Waves-Effects-On-Currents
        * Split the solution to velocity and elevation
        * Calculate the export times
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.solver_obj : The flowsolver object
        self.options    : The options of the flowsolver object
        """
        from thetis.radiation_stress import WaveEffectsOnCurrents

        # Initialise Thetis class
        ThetisTools()

        # Load dictionary with flow solver options
        dict = self.config["flowsolver 2D options"]

        print_output("- Generate Flow Solver for 2-D simulation")
        # Define Flow Solver
        self.solver_obj = solver2d.FlowSolver2d(self.mesh2d, self.bathymetry2d)

        # Calculate the start and end time of the simulation
        ThetisTools.get_start_and_end_times(self, dict)

        ## Define the flowsolver options
        self.options = self.solver_obj.options

        # Add viscosity field in solver options
        if "viscosity parameters" in self.config:
            print_output("  * Assigning horizontal eddy viscosity option")
            self.options.horizontal_viscosity = self.h_visc


        # Add weeting and drying
        if dict["use_wetting_and_drying"]:
            alpha = self.config["wetting-dry parameter"]
            self.options.wetting_and_drying_alpha = Constant(alpha)

        # Add manning field in solver options
        if "manning parameters" in self.config:
            print_output("  * Assigning manning coefficient option")
            self.options.manning_drag_coefficient = self.n_mann

        if "wind forcing" in self.config:
            print_output("  * Assigning wind stress option")
            self.options.wind_stress = self.wind_stress_2d
            # Create wind_stress_2d field to allow the user to see the
            # progression of the wind stresses
            self.solver_obj.add_new_field(
                self.wind_stress_2d,
                "wind_stress_2d",
                "Wind stress",
                "WindStress2d",
                unit="kg m-1 s-2"
            )

        # Add coriolis frequency parameter
        if "coriolis parameters" in self.config:
            print_output("  * Assigning Coriolis parameter option")
            self.options.coriolis_frequency = self.coriolis

        # Add user-defined flowsolver options
        print_output("  * Assigning user-defined solver options")
        self.options.update(dict)

        # Add user-defined timestepper options
        print_output("  * Assigning user-defined timestepper options")
        self.options.swe_timestepper_options.update(
            self.config["swe timestepper options"]
        )

        # Boundary conditions
        print_output("  * Assigning booundary conditions")
        self.solver_obj.bnd_functions['shallow_water'] = self.swe_bnd

        # Initial conditions
        print_output("  * Assigning initial conditions")
        ThetisTools.get_initial_conditions(self)
        self.solver_obj.assign_initial_conditions(
            uv=self.uv_velocities_ic,
            elev=self.elevation_ic
        )

        # Detectors
        if "detectors" in self.config:
            DetectorTools.get_detectors(self)
            # Detectors Callback
            detector_callback = DetectorsCallback(
                self.solver_obj,
                self.detector_coords,
                self.config["detectors"]["fields to export"],
                name="detectors",
                detector_names=self.detector_names
            )
            self.solver_obj.add_callback(detector_callback, 'timestep')

        # Account for Wave-Effects on currents
        # If the coupling is both ways or from SWAN-to-Thetis
        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "SWAN-to-Thetis" or
                "WEoC" in self.config
        ):
            # Create Object accounting for the Wave Effects on Currents
            self.weoc = WaveEffectsOnCurrents(
                self.mesh2d,
                self.solver_obj.options,
                self.solver_obj.fields
            )
            if "WEoC" in self.config:
                print_output("  * Add WEoC")
                self.weoc.update(
                    h_wave=self.hwave,
                    dir_wave=self.dir,
                    l_wave=self.wlen,
                    qb=self.qb,
                    c1=self.ramp_coeff
                )
            else:
                print_output("  * Add WCI")
                self.weoc.update(
                    h_wave=self.hwave_old,
                    dir_wave=self.dir_old,
                    l_wave=self.wlen_old,
                    qb=self.qb_old,
                    c1=self.ramp_coeff
                )

        # Split solution
        self.uv, self.elev = self.solver_obj.timestepper.solution.split()

        # Calculate export parameters
        ThetisTools.get_export_params(self)

        # # Setup the last things for the flowsolver
        # dump_hdf5 = self.options.export_diagnostics and not self.options.no_exports

        # TO-DO: Add the options for callbacks from iterate in solver2d
        # initial export
        self.solver_obj.t0_coupled()

        return


    def get_export_params(self):
        """
        Calculate the export times and the corresponding indices.
        Inputs:
        self

        Outputs:
        None

        Returns:
        self.export_times : A 1-D array containing the times where the specifi-
                            ed fields are exported
        self.export_ind   : 1-D array containing the corresponding indices for
                            the exportation times
        """
        # Get dictionary with flow solver parameters
        dict = self.config["flowsolver 2D options"]
        # Get export timestep
        dt = dict["simulation_export_time"]

        # Calculate the export times
        self.export_times = np.arange(self.t_start, self.t_end + dt, dt)
        # Calculate the indices of the export times
        self.export_ind = self.export_times / self.solver_obj.dt

        return


    def get_initial_conditions(self):
        """
        Get the initial conditions to be imposed in the simulation. These may
        be imported from hdf5 files or be user-defined with a constant value.
        The initial conditions are for the water elevation, the current veloci-
        ties u,v and may be for the radiation stresses (THIS NEEDS TO BE
        CHECKED)
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.elevation_ic     : A Discontinuous Galerkin function containing
                                the values of the water elevation field at the
                                start time of the simulation.
        self.u_velocity_ic    : A Discontinuous Galerkin function containing
                                the values of the u-component of the current
                                velocities field at the start time of the simu-
                                lation.
        self.v_velocity_ic    : A Discontinuous Galerkin function containing
                                the values of the u-component of the current
                                velocities field at the start time of the simu-
                                lation.
        self.uv_velocities_ic : A Discontinuous Galerkin vector containing the
                                values of the current velocities vector at the
                                start of the simulation
        """
        # Load dictionary with initial conditions parameters
        dict = self.config["initial conditions"]

        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "SWAN-to-Thetis"
        ):
            if (self.simul_kind == "ramp"):
                self.hwave_old.dat.data[:] = 0.
                self.dir_old.dat.data[:] = 0.
                self.wlen_old.dat.data[:] = 0.
                self.qb_old.dat.data[:] = 0.
            else:
                ## THIS NEEDS a CHECK to see thatwe have provided info
                dict_ = dict["import from file"]

                # Check if they are in the import_for_file dictionary
                for name in BmiThetis._ic_wave_names:
                    # If one of the names is not in the
                    if name not in dict_:
                        mes = ct.br + "Terminating Thetis: " + ct.e + \
                              "Although WCI are included in the " + \
                              "'continue' simulation, the field " + \
                              ct.bdc + name + ct.e + " is not included" + \
                              " in the " + ct.bdc + "import from file " + \
                              ct.e + "in the " + ct.bdc + "initial " + \
                              "conditions " + ct.e + "entry."
                        raise SystemExit(mes)



        # If the initial conditions are imported from hdf5 files
        if "import from file" in dict:
            # Load dictionary
            dict = dict["import from file"]

            # Load directory
            input_dir = dict["input directory"]
            # Remove input directory
            dict.pop("input directory")

            #Check if time information is in the dictionary
            if "time information" in dict:
                time_info = dict["time information"]
                dict.pop("time information")
            else:
                time_info = ""

            ## I WILL NEED TO CHECK THAT the user has specified elevation_ic
            # and uv_velocities_ic

            # Get the 'field' names of ICs
            ic_names = list(dict.keys())

            # Check that the names of the functions for water elevation and
            # current veocities are correct
            ThetisTools.check_initial_conditions_names(self, ic_names)

            ic_names = list(dict.keys())

            # Loop through the ICs to be assigned
            for name in ic_names:
                # Get the dictionary for this particular IC
                ic_dict = dict[name]

                print_output("  * Loading " + ic_dict["name"] + " field")

                ## Specify the filename
                # If time information is provided in the title
                if len(time_info)>0:
                    filename = ic_dict["name"] + "_" + str(time_info)

                # If no time information is provided in the title
                else:
                    filename = ic_dict["name"]

                # Check whether we are importing Thetis or wave ICs
                if name in self._ic_names_options:
                    # Create the function for this IC field
                    setattr(
                        self,
                        name,
                        Function(
                            getattr(self, ic_dict["function space"]),
                            name=ic_dict["name"]
                        )
                    )

                # Import the field from the h5
                IOTools.import_h5(
                    self, input_dir, filename, ic_dict["name"], getattr(self, name)
                )



        # If the initial conditions are user-defined
        else:
            # Get the 'field' names of ICs
            ic_names = list(dict.keys())

            # Check that the names of the functions for water elevation and
            # current veocities are correct
            ThetisTools.check_initial_conditions_names(self, ic_names)

            # Loop through the ICs to be assigned
            for name in ic_names:
                # Get the dictionary for this particular IC
                ic_dict = dict[name]

                print_output("  * Initialising constant " + ic_dict["name"])

                # Create the function for this IC
                setattr(
                    self,
                    name,
                    Function(
                        getattr(self,ic_dict["function space"]),
                        name=ic_dict["name"]
                    )
                )
                # Assign the value in the created function
                getattr(self, name).assign(Constant(ic_dict["value"]))

            # Convert 'scalar' velcotity functions to vector
            self.uv_velocities_ic = as_vector(
                (self.u_velocity_ic, self.v_velocity_ic)
            )


        return


    def get_start_and_end_times(self, dict):
        """
        Get the start and end time of the simulation depending on the kind of
        the simulation, i.e.
        (1) "ramp"     : It gives the current velocities, water elevation, etc.
                         to start the actual simulation we are interested. Its
                         end time should be zero.
        (2) "new"      : This simulation starts at time zero and may be preced-
                         ed by the "ramp" simulation
        (3) "continue" : This simulation starts at the time a "new" or another
                         "continue" simulation finished.
        Inputs:
        self :
        dict : A dictionary containing the flowsolver parameters

        Outputs:
        None

        Returns:

        self.t_start    : A float containing the start time of the simulation
        self.t_end      : A float containing the end time of the simulation
        self.t_duration : A float containing the duration of the simulation
        """


        # If the simulation kind is "new"
        if (self.simul_kind == "new"):
            print_output("  * Initialising new simulation")

        # If the simulation kind is ramp:
        elif (self.simul_kind == "ramp"):
            print_output("  * Initialising ramp simulation")

        elif (self.simul_kind == "ramp_cont"):
            print_output("  * Initialising ramp-continue simulation")


        # If the simulation kind is "continue":
        elif self.simul_kind == "continue":
            print_output("  * Continuing simulation")

        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "The simulation kind '" + clrtxt.bold + \
                clrtxt.darkcyan + str(self.config["simulation kind"]) + \
                clrtxt.end + "' is not valiable. Please choose one of the " + \
                "following options: " + clrtxt.bold + clrtxt.darkcyan + \
                "'ramp'" + clrtxt.end + "," + clrtxt.bold + clrtxt.darkcyan + \
                " 'new'" + clrtxt.end + "," + clrtxt.bold + clrtxt.darkcyan + \
                " 'continue'" + clrtxt.end + "."
            raise SystemExit(message)

        # Assign start time zero
        self.t_start = 0
        # Assign the user-defined end time of the simulation
        self.t_end = dict["simulation_end_time"]
        # Calculate the duration of the simulation
        self.t_duration = self.t_end - self.t_start

        return



class DetectorTools(BmiThetis):
    """
    Various methods to account for detectors in Thetis.
    This class contains the following methods:
    * create_shapefile(self, dict)
    * get_detectors(self)
    """

    def create_shapefile(self, dict):
        """
        Create a shapefile containing the detectors included in the Thetis sim-
        ulation at the specified location and with the specified EPSG by the
        user.
        Inputs:
        self :
        dict : A dictionary containing the necessary parameters to create a
               shapefile for the detectors

        Outputs:
        A shapefile containing the location of the detectors

        Returns:
        None
        """
        from fiona import collection
        from fiona.crs import from_epsg
        from pathlib import Path
        from shapely.geometry import Point, mapping

        print_output("    + Writing shapefile for detectors")

        # Create output directory if it doesn't exist
        Path(dict["output directory"]).mkdir(parents=True, exist_ok=True)

        # Determine filepath
        filepath = dict["output directory"] + dict["filename"] + ".shp"
        # If the file already exists, delete it
        Path(filepath).unlink(missing_ok=True)

        # Geometric scheme for point
        schema = {'geometry': 'Point', 'properties': {'name': 'str'}}
        # Coordinate Reference System (CRS)
        crs = from_epsg(int(dict["epsg"]))
        # Initialise the shapefile
        with collection(
                filepath, "w", "ESRI Shapefile", schema, crs=crs
        ) as output:
            # Loop through the detectors:
            for xy, name in zip(self.detector_coords, self.detector_names):
                # Create point
                point = Point(xy[0], xy[1])
                # Define properties and write to shapefile
                output.write(
                    {
                        'properties': {'name': name},
                        'geometry': mapping(point)
                    }
                )

        return


    def get_detectors(self):
        """
        From the provided csv file determine the location and the names of the
        detectors. Move any detectors that are within the maximum distance out-
        side of the domain, inside of it. If the user requests it, create a shp
        file showcasing the location of the detectors.
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.detector_coords : A list containing 2 arrays with the x- and y-co-
                               oordinates of the detectors
        self.detector_names  : A list containing the names of the detectors
        """
        # Load dictionary containing the detectors parameters
        dict = self.config["detectors"]

        print_output("  * Including detectors")

        # Get the detector names for the csvDetector names and coords from csv
        # file
        self.detector_names = np.loadtxt(
            dict["csv file"], skiprows=1, usecols=(0), dtype=str, delimiter=','
        )
        # Get the location, i.e. the coordinates, of the detectors from the csv
        # file
        self.detector_coords = np.loadtxt(
            dict["csv file"], skiprows=1, usecols=(1,2), delimiter=","
        )

        # Move any detectors that are close to the domain inside
        self.detector_coords, self.detector_names = select_and_move_detectors(
            self.mesh2d,
            self.detector_coords,
            self.detector_names,
            dict["maximum distance"]
        )

        # Check if the user has requested the creation of the shapefile
        if "create shapefile" in dict:
            DetectorTools.create_shapefile(self, dict["create shapefile"])

        return



class WCITools(BmiThetis):
    """
    Methods included in this class are:
    * check_coupling_timestep(self)
    * configure_wci(self
    * define_input_functions(self)
    * define_output_functions(self)
    * export_Thetis_coords
    * get_coupling_status(self)
    * get_coupling_timestep(self)
    * SWAN_grid(self)
    * SWAN_regular_grid(self, dict)
    """

    def check_swan_crs(self):
        """

        :return:
        """

        CoordsTools.get_mesh_coords_in_swan_crs(
            self, self.config["WCI"]["SWAN parameters"]["mesh parameters"]
        )
        # Gather the Thetis coordinates in swan CRS
        self.thetis_coords_in_swan_crs = \
            MPITools.gather_provided_coords_across_ranks(
                self, self.swan_coords
            )

        return



    def check_coupling_timestep(self):
        """
        Check that the coupling and simulation timestep are either multiples or
        factors of one another. If they aren't, terminate Thetis
        Inputs:
        self :

        Outputs:
        None

        Returns:
        None
        """
        # Get the simulation timestep from the solver options
        sim_dt = self.config["flowsolver 2D options"]["timestep"]

        # If the simulation timestep is bigger than the coupling timestep
        if sim_dt >= self.coupl_dt:
            # If the simulation timestep is not a multiple of the coupling ti-
            # mestep, terminate Thetis
            if sim_dt%self.coupl_dt != 0:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The simulation timestep " + clrtxt.bold + \
                    clrtxt.darkcyan + str(sim_dt) + clrtxt.end + " s is " + \
                    "not a multiple of the coupling timestep" + clrtxt.bold + \
                    clrtxt.darkcyan + " " + str(self.coupl_dt) + clrtxt.end + \
                    " s."
                raise SystemExit(message)

        # If the coupling timestep is bigger than the simulation timestep
        else:
            # If the coupling timestep is not a multiple of the simulation ti-mestpe, terminate Thetis
            if self.coupl_dt%sim_dt != 0:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The coupling timestep " + clrtxt.bold + \
                    clrtxt.darkcyan + str(self.coupl_dt) + clrtxt.end + \
                    " s is not a multiple of the simulation timestep " + \
                    clrtxt.bold + clrtxt.darkcyan + str(sim_dt) + \
                    clrtxt.end + " s."
                raise SystemExit(message)

        return


    def configure_wci(self):
        """
        Configure Wave-Current Interactions (WCI) between Thetis and SWAN. To
        do so:
        * determine the coupling status: 2-way, Thetis-to-SWAN or SWAN-to-Thetis
        * determine the coupling timestep
        * if the simulation is run in parallel,
          - make the mesh coordinates available across all ranks
          - determine the starting and ending indices in the global array for
            the rank arrays
        * if 2-way ir SWAN-to-Thetis
          - export a txt with the Thetis coords for SWAN
          - define the input functions to put the variables from SWAN
        * if the coupling is 2-way or from Thetis-to-SWAN
          - calculate the coordinates of the SWAN grid
          - define the functions for the Thetis fields that will be exported to
            SWAN
        Inputs:
        self :

        Outputs:
        None

        Returns:
        None
        """
        from scipy.interpolate import NearestNDInterpolator, griddata

        # # Get the coupling status between the models
        # WCITools.get_coupling_status(self)

        # Get the coupling timestep
        WCITools.get_coupling_timestep(self)

        # Check that the coupling timestep is a multiple of the simulation ti-
        # mestep or that the simulation timestep is a multiple of the coupling
        # timestep
        WCITools.check_coupling_timestep(self)


        # Check the coordinates system
        WCITools.check_swan_crs(self)



        # If the coupling is both ways or from SWAN-to-Thetis
        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "SWAN-to-Thetis"
        ):
            # Export a txt with the Thetis coordinates for SWAN
            WCITools.export_Thetis_coords(self)

            # Define functions for input variables from SWAN
            WCITools.define_input_functions(self)

            # Check if the timestep of Thetis is smaller than the the coupling
            # one
            dt_th = self.config["flowsolver 2D options"]["timestep"]
            if (self.coupl_dt - dt_th) > 10 ** (-8):
                # We will need to interpolate the wave fields in time
                self.wave_char_interp = True
                self.no_wave_interp = int(self.coupl_dt / dt_th)

                # Initialise funcstionfor the differences


            else:
                self.wave_char_interp = False

        # If the coupling is both ways or from Thetis-to-SWAN
        if (
                self.coupl_stat == "2-way" or
                self.coupl_stat == "Thetis-to-SWAN"
        ):
            # Construct the SWAN grid
            WCITools.SWAN_grid(self)

            # Get the elevationa dn velocity in CG
            self.elev_CG = Function(self.P1_2d, name='CG_elevation_2d')
            self.uv_CG = Function(self.P1v_2d, name="CG_velocities_2d")

            # Get the minimum depth
            self.swan_min_dep = \
                self.config["WCI"]["SWAN parameters"]["threshold depth [m]"]

            # Define functions for Thetis output variables
            WCITools.define_output_functions(self)

            # Calculate the bathymetry at SWAN grid points
            bath_data = self.comm.allgather(self.bathymetry2d.dat.data)
            bath_data = np.concatenate(bath_data)

            if self.mpi_stat:
                self.bath_data_at_swan = griddata(
                    self.xy_gathered,
                    bath_data,
                    (
                        self.swan_coords_in_mesh_crs[:, 0],
                        self.swan_coords_in_mesh_crs[:, 1]
                    ),
                    method='linear',
                    fill_value=0
                )
            else:
                self.bath_data_at_swan = griddata(
                    self.xy,
                    bath_data,
                    (
                        self.swan_coords_in_mesh_crs[:, 0],
                        self.swan_coords_in_mesh_crs[:, 1]
                    ),
                    method='linear',
                    fill_value=0
                )

            # Find the indices where we have positive bathymetry less than the
            # threshold value
            ind = np.where(
                ( self.bath_data_at_swan[:] > 0.00 ) &
                ( self.bath_data_at_swan[:] < self.swan_min_dep)
            )[0]
            self.bath_data_at_swan[ind] = self.swan_min_dep

            # Calculate the indices where we are smaller than zero, i.e. we
            # have land.We decided to impose a cut-off limit of 0.02 m, i.e.
            # for bathymetry that is less than 0.05 m, we consider it land
            self.ind_land = np.where(
                self.bath_data_at_swan[:] < 0.05
            )[0]

        return


    def define_input_functions(self):
        """
        Create functions for any variable that can be imported from SWAN. The
        names of the functions are taken from the corresponding output name the
        variable has in SWAN (to avoid confusion)
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.hwave   : A function for the wave height
        self.dir  : A function for the mean wave direction
        self.wlen : A function for the mean wavelength
        self.qb   : A function for the percentage of wave breaking
        """
        # Define array for the wave height
        self.hwave = Function(self.P1_2d, name="wave_height")
        self.hwave_old = Function(self.P1_2d, name="old_wave_height")
        self.hwave_new = Function(self.P1_2d, name="new_wave_height")
        self.hwave_diff = Function(
            self.P1_2d, name="wave_height_differences"
        )

        # Define array for the mean wave direction
        self.dir = Function(self.P1_2d, name="mean_wave_direction")
        self.dir_old = Function(self.P1_2d, name="old_mean_wave_direction")
        self.dir_new = Function(self.P1_2d, name="new_mean_wave_direction")
        self.dir_diff = Function(
            self.P1_2d, name="mean_wave_direction_differences"
        )

        # Define array for the mean wavelength
        self.wlen = Function(self.P1_2d, name="mean_wavelength")
        self.wlen_old = Function(self.P1_2d, name="old_mean_wavelength")
        self.wlen_new = Function(self.P1_2d, name="new_mean_wavelength")
        self.wlen_diff = Function(
            self.P1_2d, name="mean_wavelength_differences"
        )

        # Define array for the percentage of wave-braking
        self.qb = Function(self.P1_2d, name="wave_breaking_percentage")
        self.qb_old = Function(self.P1_2d, name="old_wave_breaking_percentage")
        self.qb_new = Function(self.P1_2d, name="new_wave_breaking_percentage")
        self.qb_diff = Function(
            self.P1_2d, name="wave_breaking_percentage_differences"
        )

        return


    def define_output_functions(self):
        """
        Create functions for any Thetis variable/field that can be exported and
        imported into SWAN. The names of the variables are the names of the
        corresponding variables in the SWAN code
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.wlevl : A function for the water level
        self.uxb   : A function for the x-component of the current field
        self.uyb   : A function for the y-component of the current field
        """
        # Define function for water level
        self.wlevl = np.empty(np.shape(self.xy_swan)[0])
        # Define function for the x-component of the current field
        self.uxb = np.empty(np.shape(self.xy_swan)[0])
        # Define function for the y-component of the current field
        self.uyb = np.empty(np.shape(self.xy_swan)[0])

        return


    def export_Thetis_coords(self):
        """
        Export the Thetis mesh coordinates in a temporary directory to be read
        by SWAN.
        Inputs:
        self :

        Outputs:
        A txt file containing the mesh coordinates in a temporary directory

        Returns:
        None
        """
        from pathlib import Path

        # Specify temporary output directory
        output_dir = './temp/'
        # Create directory if it doesn't exist
        Path(output_dir).mkdir(parents=True, exist_ok=True)


        # If the simulation is run in parallel
        if self.mpi_stat:
            # For the 0th rank only
            if self.mpi_rank==0:
                # Save coords to txt
                np.savetxt(
                    output_dir + 'coords.txt',
                    self.thetis_coords_in_swan_crs,
                    fmt="%s",
                    delimiter=" ",
                    newline='\n'
                )

            self.comm.Barrier()

        # If the simulation is run serially
        else:
            # Save coords to txt
            np.savetxt(
                output_dir + 'coords.txt',
                self.thetis_coords_in_swan_crs,
                fmt="%s",
                delimiter=" ",
                newline='\n'
            )

        return


    def get_coupling_status(self):
        """
        Get the coupling status between Thetis and SWAN. Acceptable options are
        '2-way', 'SWAN-to-Thetis', and 'Thetis-to-SWAN'.
        Inputs:
        self

        Outputs:
        None

        Returns:
        self.coupl_stat : A string containing the coupling status between Thet-
                          is and SWAN
        """
        # Load the dictionary with the WCI parameters
        dict = self.config["WCI"]

        # If the models are fully coupled
        if dict["coupling"] == "2-way":
            self.coupl_stat = "2-way"
        elif dict["coupling"] == "SWAN-to-Thetis":
            self.coupl_stat = "SWAN-to-Thetis"
        elif dict["coupling"] == "Thetis-to-SWAN":
            self.coupl_stat = "Thetis-to-SWAN"
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      clrtxt.end + "The coupling option " + clrtxt.bold + \
                      clrtxt.darkcyan + "'" + str(dict["coupling"]) + "' " + \
                      clrtxt.end + "is not a valid option. Use one of the " + \
                      "following: " + clrtxt.bold + clrtxt.darkcyan + "'2-way'" + \
                      clrtxt.end + ", " + clrtxt.bold + clrtxt.darkcyan + \
                      "'SWAN-to-Thetis'" + clrtxt.end + ", " + clrtxt.bold + \
                      clrtxt.darkcyan + "'Thetis-to-SWAN'" + clrtxt.end
            raise SystemExit(message)

        return


    def get_coupling_timestep(self):
        """
        Get the coupling timestep in seconds
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.coupl_dt : The coupling timestep between the models [s]
        """
        # Load the dictionary with the WCI parameters
        dict = self.config["WCI"]

        # If the timestep hasn't been included in the WCI parameters, terminate
        # Thetis
        if "timestep [s]" not in dict:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "No coupling timestep has been defined. Add " + \
                "the entry " + clrtxt.bold + clrtxt.darkcyan + "'timestep" + \
                " [s]' " + clrtxt.end + "under the entry " + clrtxt.bold + \
                clrtxt.darkcyan + "'WCI'" + clrtxt.end
            raise SystemExit(message)

        self.coupl_dt = dict["timestep [s]"]

        return


    def SWAN_grid(self):
        """
        Depending on the SWAN grid type, calculate the coordinates of SWAN by
        calling the appropriate method.
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.xy_swan : 2-D array containing the coordinates of SWAN mesh
        """
        # Load WCI parameters from json file
        dict = self.config["WCI"]

        # If there is no 'SWAN parameters' entry, terminate Thetis
        if "SWAN parameters" not in dict:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "No entry " + clrtxt.bold + clrtxt.darkcyan + \
                "'SWAN parameters' " + clrtxt.end + "exists under the " + \
                "entry " + clrtxt.bold + clrtxt.darkcyan + "'WCI'" + clrtxt.end
            raise SystemExit(message)

        # If there is no "mesh parameters" entry, terminate Thetis
        if "mesh parameters" not in dict["SWAN parameters"]:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "No entry " + clrtxt.bold + clrtxt.darkcyan + \
                "'mesh parameters' " + clrtxt.end + "exists under the " + \
                "entry " + clrtxt.bold + clrtxt.darkcyan + "'SWAN " + \
                "parameters' " + clrtxt.end
            raise SystemExit(message)

        # Load SWAN's mesh parameters
        dict = dict["SWAN parameters"]["mesh parameters"]

        # If the SWAN mesh is regular
        if dict["type"] == "Regular":
            # Calculate the coordinates
            WCITools.SWAN_regular_grid(self, dict)

        # Else if the SWAN mesh is not regular
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                "When coupling with SWAN, the SWAN mesh type " + \
                str(dict["type"]) + " has not been programmed yet" + clrtxt.end
            raise SystemExit(message)

        CoordsTools.get_swan_coords_in_mesh_crs(self, dict)

        return


    def SWAN_regular_grid(self, dict):
        """
        Calculate the orthogonal structured SWAN mesh coordinates according to
        the parameters specified in the json file.
        Inputs:
        self :
        dict : A dictionary containing the mesh parameters of SWAN

        Outputs:
        None

        Returns:
        self.xy_swan : 2-D array containing the coordinates of SWAN mesh
        """
        # Calculate the values of the x- and y-coordinates
        x = np.linspace(
            dict["xpc"], dict["xpc"] + dict["xlenc"], dict["mxc"] + 1
        )
        y = np.linspace(
            dict["ypc"], dict["ypc"] + dict["ylenc"], dict["myc"] + 1
        )

        # Grid coordinates
        xv, yv = np.meshgrid(x, y)

        # Flatten
        x = xv.flatten()
        y = yv.flatten()

        # Stack them in an array
        self.xy_swan = np.stack((x, y), axis=-1)

        return









    def import_fields_from_SWAN(self):
        """
        """
        if self.config["WEC"]["Thetis_calculations"]:
            print_output("     + Calculate WEC in Thetis")

    def calculate_output_fields(self):
        """
        """
        from scipy.interpolate import NearestNDInterpolator, griddata

        # Convert elevation to CG
        # self.elev_CG.project(self.elev)
        self.elev_CG.interpolate(self.elev)


        # Gather and broadcast the values across all ranks
        elev_CG_data = self.comm.allgather(self.elev_CG.dat.data)
        elev_CG_data = np.concatenate(elev_CG_data)

        # Create interpolator
        # intp = NearestNDInterpolator(self.xy_gathered, elev_CG_data)
        # If the run is serial
        if self.mpi_stat:
            data = griddata(
                self.xy_gathered,
                elev_CG_data,
                (
                    self.swan_coords_in_mesh_crs[:, 0],
                    self.swan_coords_in_mesh_crs[:, 1]
                ),
                method='linear',
                fill_value=0
            )
        else:
            data = griddata(
                self.xy,
                elev_CG_data,
                (
                    self.swan_coords_in_mesh_crs[:, 0],
                    self.swan_coords_in_mesh_crs[:, 1]
                ),
                method='linear',
                fill_value=0
            )

        self.wlevl = data * self.ramp_coeff

        # Velocities
        # self.uv_CG.project(self.uv)
        self.uv_CG.interpolate(self.uv)

        # Gather and broadcast the values across all ranks
        uv_CG_data = self.comm.allgather(self.uv_CG.dat.data)
        uv_CG_data = np.concatenate(uv_CG_data)

        # Interpolate
        if self.mpi_stat:
            data = griddata(
                self.xy_gathered,
                uv_CG_data,
                (
                    self.swan_coords_in_mesh_crs[:, 0],
                    self.swan_coords_in_mesh_crs[:, 1]
                ),
                method='linear',
                fill_value=0
            )
        else:
            data = griddata(
                self.xy,
                uv_CG_data,
                (
                    self.swan_coords_in_mesh_crs[:, 0],
                    self.swan_coords_in_mesh_crs[:, 1]
                ),
                method='linear',
                fill_value=0
            )

        self.uxb = data[:, 0] * self.ramp_coeff
        self.uyb = data[:, 1] * self.ramp_coeff


        ## Zero the values of eta and uv where we have 'land' originally
        if np.shape(self.ind_land)[0] > 0:
            self.uxb[self.ind_land] = 0.0
            self.uyb[self.ind_land] = 0.0
            self.wlevl[self.ind_land] = 0.0

        ## Calculate any indices where the water dries
        # Calculate the water depth
        water_depth = self.wlevl + self.bath_data_at_swan
        # FInd the indices where the water depth is smaller than 0.05 m and it
        # isn't land
        ind_dry = np.where(
            (water_depth < 0.05) & (self.bath_data_at_swan >= 0.05)
        )[0]
        # Check if we have dried land
        if np.shape(ind_dry)[0]  != 0:
            # print("We have dry land")
            self.wlevl[ind_dry] = np.fmin(
                np.zeros(np.shape(ind_dry)), self.wlevl[ind_dry] + 0.05
            )
            # # if we zeored out the water elevaion
            # ind_dry = np.where(self.wlevl[ind_dry] == 0 )[0]
            # self.uxb[ind_dry] = 0.0
            # self.uyb[ind_dry] = 0.0

        #     print(self.wlevl[ind_dry])
        #     print("----ffffffffffffffffffffffffffffffffffff")
        #
        # print(f"u: [{np.amin(self.uxb):.3f}, {np.amax(self.uxb):.3f}] m/s")
        # print(f"v: [{np.amin(self.uyb):.3f}, {np.amax(self.uyb):.3f}] m/s")
        # print(f"eta: [{np.amin(self.wlevl):.3f}, {np.amax(self.wlevl):.3f}] m")

        # print(self.ind_land)
        # print(" ")
        # print(np.shape(ind_dry))
        # print(self.bath_data_at_swan[ind_dry])
        # print(self.wlevl[ind_dry])


        # print(" ")
        # print(ind_dry)
        # exit()

        # Keep the velocity, but minimise the water elevation to
        # if np.shape(ind_dry)[0] > 0:
        #     self.uxb[ind_dry] = 0.0
        #     self.uyb[ind_dry] = 0.0
        #     self.wlevl[ind_dry] = 0.0
            # # Find the indices where the water becomes too shallow
            # data = self.wlevl + self.bath_data_at_swan
            # ind_dry = np.where(data < 0)[0]
            # print(ind_dry)
            # print(ind_land)
            #
            # exit()

        return

    def calculate_fields_for_SWAN(self, f_name, grid_name):
        """
        Calculate values at mesh indices (either Thetis or SWAN's mesh) for
        either the water elevation (elev) or currents (uv) to be imported into
        SWAN

        Inputs:
        f_name    : The name of the function to be evaluated
        grid_name : The name of the grid to be evaluated upon.
                    Two options available: Thetis or SWAN

        Outputs:
        data : The values of the function at the specified grid. Flattened array
        """
        from scipy.interpolate import NearestNDInterpolator

        # Initialize lists for interpolator
        xy_intp = []
        z_intp = []
        # Get coordinates and values for interpolator
        for index in range(np.shape(self.xy_gathered)[0]):
            value = getattr(self, f_name).at(self.xy_gathered[index],
                                             dont_raise=True)
            if value is not None:
                xy_intp.append([
                    self.xy_gathered[index, 0],
                    self.xy_gathered[index, 1]
                ])
                z_intp.append(value)
        z_intp = np.concatenate(z_intp, axis=None)
        # Create interpolator
        intp = NearestNDInterpolator(xy_intp, z_intp)

        if grid_name == self._name:
            data = intp(self.xy_gathered)
        elif grid_name == "SWAN":
            data = intp(self.xy_swan)
        else:
            print_output(f"No grid '{grid_name}' exists." \
                         " Terminating simulation."
                         )
            exit()

        if f_name == "uv":
            data = np.transpose(data)

        return data

    def set_input_field(self, name, data):
        """
        """
        if name == 'sea_surface_water_wave__height':
            if self.mpi_stat:
                self.hwave_new.dat.data[:] = data[
                    self.mpi_start_ind[self.mpi_rank]:
                            self.mpi_end_ind[self.mpi_rank]
                    ]
            else:
                self.hwave_new.dat.data[:] = data



        elif name == 'sea_surface_water_wave__direction':

            if self.mpi_stat:
                self.dir_new.dat.data[:] = data[
                                               self.mpi_start_ind[self.mpi_rank]:
                                               self.mpi_end_ind[self.mpi_rank]
                                            ]
            else:
                self.dir_new.dat.data[:] = data

        elif name == 'sea_surface_water_wave__wavelength':

            if self.mpi_stat:
                self.wlen_new.dat.data[:] = data[
                                                self.mpi_start_ind[self.mpi_rank]:
                                                self.mpi_end_ind[self.mpi_rank]
                                                ]
            else:
                self.wlen_new.dat.data[:] = data


        elif name == 'sea_surface_water_wave__breaking_fraction':

            if self.mpi_stat:
                self.qb_new.dat.data[:] = data[
                                            self.mpi_start_ind[self.mpi_rank]:
                                            self.mpi_end_ind[self.mpi_rank]
                                            ]
            else:
                self.qb_new.dat.data[:] = data


        return




    def get_output_field(self, name):
        if name == 'sea_water_surface__elevation':
            src = self.wlevl
        elif name == 'sea_water_flow__x_component_of_velocity':
            src = self.uxb
        elif name == 'sea_water_flow__y_component_of_velocity':
            src = self.uyb

        return src

    # def get_indices_per_rank(self):
    #     """
    #     Find the starting index as well as the length of the data per rank
    #     """
    #     import numpy as np
    #     # Number of elements per rank
    #     self.counts = self.comm.allgather(
    #         np.shape(
    #             self.xy
    #         )[0]
    #     )
    #     # Starting displacement, i.e. starting indices for the global array
    #     self.s_index = [0] + list(np.cumsum(self.counts)[0:-1])
    #     # Convert to numpy arrays
    #     self.counts = np.asarray(self.counts)
    #     self.s_index = np.asarray(self.s_index)
    #     # Ending indices per rank with respect to the global array
    #     self.e_index = self.s_index + self.counts

    # def thetis_coords_across_ranks(self):
    #     """
    #     """
    #
    #     # # Separate coordinates
    #     # xy =
    #
    #     self.xy_gathered = np.concatenate(xy)
    #
    #
    #     # self.xy_gathered = np.around(self.xy_gathered, 4)

    # def swan_structured_coords(self):
    #     """
    #     """
    #     # x- and y- coords
    #     x = np.linspace(
    #         self.config["WEC"]["mesh parameters"]["xpc"],
    #         self.config["WEC"]["mesh parameters"]["xpc"] + \
    #         self.config["WEC"]["mesh parameters"]["xlenc"],
    #         self.config["WEC"]["mesh parameters"]["mxc"] + 1
    #     )
    #     y = np.linspace(
    #         self.config["WEC"]["mesh parameters"]["ypc"],
    #         self.config["WEC"]["mesh parameters"]["ypc"] + \
    #         self.config["WEC"]["mesh parameters"]["ylenc"],
    #         self.config["WEC"]["mesh parameters"]["myc"] + 1
    #     )
    #
    #     # Grid coordinates
    #     xv, yv = np.meshgrid(x, y)
    #
    #     # Flatten
    #     x = xv.flatten()
    #     y = yv.flatten()
    #
    #     # Stack them in an array
    #     self.xy_swan = np.stack((x, y), axis=-1)






class TideTools(BmiThetis):
    """
    Various methods for tide forcing.
    This class contains the following methods:
    * configure_tide_forcing(self)
    * configure_tide_interp(self)
    * create_tide_function_vs_time(self, dict)
    * determine_tide_start_date(self, dict)
    * set_tidal_field(self, t)
    * update_tide_time_function(self, t, item)
    """

    def __init__(self):
        BmiThetis._tide_options = [
            "constituent forcing",
            "TPXO forcing",
            "AMCG forcing"
        ]


    def configure_tide_forcing(self):
        """
        Depending on how the tide constituents are provided, the appropriate
        method for configuring the tides is selected
        Inputs:
        self

        Outputs:
        NOne

        Returns:
        None
        """
        # Initiliase TideClass
        TideTools()

        # Load dictionary regarding the tide forcing
        dict = self.config["tide forcing"]

        if list(dict.keys())[0] not in self._tide_options:
            mes = ct.br + "Termintaing Thetis: " + ct.e + "The option " + \
                  ct.bdc + str(list(dict.keys())[0]) + ct.e + " is not " + \
                  "a viable option for tide forcing. Please choose one of " + \
                  "the following:" + ct.bdc + str(self._tide_options)[1:-1] + \
                  ct.e
            raise SystemExit(mes)

        # If we provide the constituents through a csv and not from the TPXO
        if "constituent forcing" in dict:
            TideTools.create_tide_function_vs_time(
                self, dict["constituent forcing"]
            )

        # If the constituents are from the TPXO
        elif "TPXO forcing" in dict:
            TideTools.configure_tpxo_tide_interp(self, dict["TPXO forcing"])

        # If the constituents are in AMCG format
        elif "AMCG forcing" in dict:
            TideTools.configure_amcg_tide_interp(self, dict["AMCG forcing"])

        return


    def configure_amcg_tide_interp(self, dict):
        """

        :param dict:
        :return:
        """
        from uptide import Tides, tidal_netcdf

        print_output("  * Configuring tide AMCG interpolator")

        # Determine harmonic constituents to be included in the tidal computa-
        # tions
        tide = Tides(dict["constituents"])

        # Get start datetime information
        start_date = TideTools.determine_tide_start_date(self, dict)

        # Set t=0 for the start date of the tides
        tide.set_initial_time(start_date)

        # # Create a tuple containing the limits of the interplator
        # coord_lim = (
        #     tuple(dict["longitude limits [o]"]),
        #     tuple(dict["latitude limits [o]"])
        # )

        # Create tide constitudents interpolator
        self.tide_interp = tidal_netcdf.AMCGTidalInterpolator(
            tide, dict["nc file"]
        )

        if "boundary_ids" in dict:
            # Get the TPXO coords:
            CoordsTools.get_mesh_coords_in_TPXO_crs(self, dict["boundary_ids"])
        else:
            CoordsTools.get_mesh_coords_in_TPXO_crs(self)

        return


    def configure_tpxo_tide_interp(self, dict):
        """
        Create a tide interpolator with t=0 the start date of the tides defined
        in the json file with as many tide constituents as the user specifies
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.tide_interp : An interpolator of tides for the specified harmonic
                           constituents
        """
        from uptide import Tides, tidal_netcdf


        print_output("  * Configuring tide TPXO interpolator")

        # Determine harmonic constituents to be included in the tidal computa-
        # tions
        tide = Tides(dict["constituents"])

        # Get start datetime information
        start_date = TideTools.determine_tide_start_date(self, dict)

        # Set t=0 for the start date of the tides
        tide.set_initial_time(start_date)

        # Convert longitude coordinates from [-180,180] to [0,360] since the
        # TPXO has longitude with range [0,360]
        long_lim = dict["longitude limits [o]"]
        # long_lim = [i+180 for i in long_lim]
        # Create a tuple containing the limits of theinterpolator
        coord_lim = (
            tuple(long_lim),
            tuple(dict["latitude limits [o]"])
        )

        # Create tide constituents interpolator
        self.tide_interp = tidal_netcdf.TPXOTidalInterpolator(
            tide,
            dict["grid file"],
            dict["data file"],
            ranges = coord_lim
        )

        if "boundary_ids" in dict:
            # Get the TPXO coords:
            CoordsTools.get_mesh_coords_in_TPXO_crs(
                self, dict["boundary_ids"], TPXO_forcing=True
            )
        else:
            CoordsTools.get_mesh_coords_in_TPXO_crs(self, TPXO_forcing=True)


        return


    def configure_amcg_uv_interp(self, dict):
        """

        """
        from uptide import Tides, tidal_netcdf

        print_output("  * Configuring uv tidal AMCG interpolator")

        # Determine harmonic constituents to be included in the tidal computa-
        # tions
        tide = Tides(dict["constituents"])

        # Get start datetime information
        start_date = TideTools.determine_tide_start_date(self, dict)

        # Set t=0 for the start date of the tides
        tide.set_initial_time(start_date)

        # # Create a tuple containing the limits of the interplator
        # coord_lim = (
        #     tuple(dict["longitude limits [o]"]),
        #     tuple(dict["latitude limits [o]"])
        # )

        # Create u-velocity constituents interpolator
        self.uv_x_int = tidal_netcdf.AMCGTidalInterpolator(
            tide, dict["u nc file"]
        )
        self.uv_y_int = tidal_netcdf.AMCGTidalInterpolator(
            tide, dict["v nc file"]
        )

        if "boundary_ids" in dict:
            # Get the TPXO coords:
            CoordsTools.get_mesh_coords_in_TPXO_crs(self, dict["boundary_ids"])
        else:
            CoordsTools.get_mesh_coords_in_TPXO_crs(self)

        # Check that we have uv_interpolator in the BC
        ForcingTools.check_vector_uv_update_forcings(self)



    def create_tide_function_vs_time(self, dict):
        """
        Initialise the tide function by setting the tide constituents and the
        start data. Sort out the appropriate amplitudes and phases based on the
        constituents to be used.
        Inputs:
        self :
        dict : Dictionary containing the parameters for the constituent forcing

        Outputs:
        None

        Returns:
        self.tide       : The tidal computational object
        self.tide_ampl  : 1-D array containing the amplitudes of the harmonic
                          constituents [m]
        self.tide_phase : 1-D array containing the phases of the harmonic cons-
                          tituents [rad]
        """
        import uptide

        # Load all constituents information
        names = IOTools.read_data_from_csv(
            self, dict["csv file"], 0, str, rows=1
        )
        # Load the constituents amplitude [m]
        ampl = IOTools.read_data_from_csv(
            self, dict["csv file"], 1, np.float64, rows=1
        )
        # Load the constituents phase in degrees
        phase = IOTools.read_data_from_csv(
            self, dict["csv file"],2, np.float64, rows=1
        )

        # Get the indices of the constituentd we want
        ind = np.nonzero(np.in1d(names, dict["constituents"]))[0]
        # Get the corresponding amplitudes
        self.tide_ampl = ampl[ind]
        # Get the corresponding phases
        self.tide_phase = phase[ind]
        # Convert degrees to radians
        self.tide_phase = np.radians(phase)

        # Determine harmonic constituents to be included in the tidal computa-
        # tions
        self.tide = uptide.Tides(names[ind])

        # Get start datetime information
        start_date = TideTools.determine_tide_start_date(self, dict)
        # Set t=0 for the start date of the tides
        self.tide.set_initial_time(start_date)

        return


    def determine_tide_start_date(self, dict):
        """
        Create a datetime variable with timezone information depending on the
        information on the json file, which will determine the appropriate met-
        hod to be called.
        Inputs:
        self :
        dict : A dictionary containing the parameters for TPXO forcing

        Outputs:
        None

        Returns:
        start_date : The datetime of the tides start date
        """
        # Check to see if a timezone has been specified, else terminate Thetis
        if "timezone" not in dict:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "No timezone has been specified for the tides" + \
                " configuration. Add the entry " + clrtxt.bold + \
                clrtxt.darkcyan + "'timezone' " + clrtxt.end + "under the " + \
                "entry " + clrtxt.bold + clrtxt.darkcyan + "'TPXO " + \
                "forcing' " + clrtxt.end
            raise SystemExit(message)

        # If the timezone is UTC
        if dict["timezone"] == "UTC":
            # Calculate the UTM datetime
            start_date = TimeTools.create_UTC_date(self, dict)

        # Else, terminate Thetis
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                "Other timezones apart from UTC have not been programmed " + \
                "yet for the tides start date" + clrtxt.end
            raise SystemExit(message)

        return start_date


    def set_tidal_field(self, t, item):
        """
        Set the tidal field based on the TPXO netcdf
        Inputs:
        self :
        t    : The time to calculate the forcing
        item : The entry (key, value) in the 'update forcings" entry

        Outputs:
        None

        Returns:
        The update tidal function
        """
        from uptide import netcdf_reader

        # Set the timestamp for the tide interpolator
        self.tide_interp.set_time(t + self.t_start)

        # Get the function representing the TPXO forcing
        func = getattr(self, item[0])
        # Set up the elevation function as editable field
        data = func.dat.data[:]
        data[:] = 0.

        # print(self.mpi_rank, np.shape(data))
        # print(data)
        #
        # # print(self.mpi_rank, self.bc_ind, no_owned)
        # print(" ")
        # exit()

        # Loop through the mesh-TPXO coordinates
        for i in range(np.shape(self.tpxo_coords)[0]):
            # If the point is in the sea
            try:
                data[self.bc_ind[i]] = self.tide_interp.get_val(
                    (
                        self.tpxo_coords[i,1],
                        self.tpxo_coords[i,0]
                    )
                )*self.ramp_coeff
            except netcdf_reader.CoordinateError:
                data[self.bc_ind[i]] = 0.

        #print(f"[{min(data):.3f}, {max(data):.3f}]")
        # Update the values to the function
        func.dat.data[:] = data

        return



    def set_u_tidal_field(self, t, item):
        """
        Set the tidal field based on the TPXO netcdf
        Inputs:
        self :
        t    : The time to calculate the forcing
        item : The entry (key, value) in the 'update forcings" entry

        Outputs:
        None

        Returns:
        The update tidal function
        """
        from uptide import netcdf_reader

        # Set the timestamp for the tide interpolator
        self.uv_x_int.set_time(t + self.t_start)

        # Get the function representing the TPXO forcing
        func = getattr(self, item)
        # Set up the elevation function as editable field
        data = func.dat.data[:]
        data[:,:] = 0.


        # Loop through the mesh-TPXO coordinates
        for i in range(np.shape(self.tpxo_coords)[0]):

            # If the point is in the sea
            try:
                data[self.bc_ind[i],0] = self.uv_x_int.get_val(
                    (
                        self.tpxo_coords[i,1],
                        self.tpxo_coords[i,0]
                    )
                )*self.ramp_coeff
            except netcdf_reader.CoordinateError:
                data[self.bc_ind[i], 0] = 0.

        #print(f"[{min(data):.3f}, {max(data):.3f}]")
        # Update the values to the function
        func.dat.data[:] = data

        return

    def set_v_tidal_field(self, t, item):
        """
        Set the tidal field based on the TPXO netcdf
        Inputs:
        self :
        t    : The time to calculate the forcing
        item : The entry (key, value) in the 'update forcings" entry

        Outputs:
        None

        Returns:
        The update tidal function
        """
        from uptide import netcdf_reader

        # Set the timestamp for the tide interpolator
        self.uv_y_int.set_time(t + self.t_start)

        # Get the function representing the TPXO forcing
        func = getattr(self, item)
        # Set up the elevation function as editable field
        data = func.dat.data
        data[:,:] = 0.

        # Loop through the mesh-TPXO coordinates
        for i in range(np.shape(self.tpxo_coords)[0]):
            # If the point is in the sea
            try:
                data[i, 1] = self.uv_y_int.get_val(
                    (
                        self.tpxo_coords[i, 1],
                        self.tpxo_coords[i, 0]
                    )
                ) * self.ramp_coeff
            except netcdf_reader.CoordinateError:
                data[i, 1] = 0.

        # print(f"[{min(data):.3f}, {max(data):.3f}]")
        # Update the values to the function
        func.dat.data[:] = data

        return



    def update_tide_time_function(self, t, item):
        """
        Calculate the water elevation based on the given amplitudes and phases
        for time t and assign this value to the appropriate scalar function
        Inputs:
        t    : The time of the simulation
        itme : he entry (key, value) in the 'update forcings" entry

        Outputs:
        None

        Returns:
        None
        """

        # Calculate the water elevation
        value = self.tide.from_amplitude_phase(
            self.tide_ampl, self.tide_phase, t
        ) * self.ramp_coeff

        # Assign the value to the function
        getattr(self, item[0]).assign(Constant(value))

        return






class MPITools(BmiThetis):
    """
    A class to aid with any actions dealing with parallel run
    The methods contained in this class are:
    * gather_coords_across_ranks(self)
    * get_indices_per_rank(self)
    * get_mpi_status(self)

    """
    def gather_coords_across_ranks(self):
        """
        Gather the mesh coordinates and have them available across all ranks
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.xy_gathered : 2-D array containing all the mesh coordinates in a
                           parallel run
        """
        # Gather all the coordinates across ranks
        xy = self.comm.allgather(self.xy)
        # Flatten to a 2-D array
        self.xy_gathered = np.concatenate(xy)

        return

    def gather_provided_coords_across_ranks(self, coords):
        """

        :param coords:
        :return:
        """
        xy = self.comm.allgather(coords)
        xy_gathered = np.concatenate(xy)

        return xy_gathered



    def get_indices_per_rank(self):
        """
        If the model is run in parallel, calculate the number of elements each
        rank has, their starting index in the gathered array, as well as their
        ending index + 1 (perfect for array slicing) in the gathered array
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.mpi_rank_len  : A 1-D array containing the amount of elements per
                             rank
        self.mpi_start_ind : A 1-D array containing the starting index of the
                             rank's array elements in the gathered array
        self.mpi_end_ind   : A 1-D array containing the final index + 1 of the
                             rank's array elements in the gathered array, which
                             will be used for array slicing
        """
        # Number of elements per rank
        self.mpi_rank_len = self.comm.allgather(np.shape(self.xy)[0])

        # Starting index per rank for the global array
        self.mpi_start_ind = [0] + list(np.cumsum(self.mpi_rank_len)[0:-1])

        # Convert list to arrays
        self.mpi_rank_len = np.asarray(self.mpi_rank_len)
        self.mpi_start_ind = np.asarray(self.mpi_start_ind)

        # Final index per rank for the global array for array/list slicing
        self.mpi_end_ind = self.mpi_start_ind + self.mpi_rank_len

        return


    def get_mpi_status(self):
        """
        Get basic MPI information, like the MPI communicator, the number or mpi
        processes, the ID of the rank. Determine if we run in parallel or not.
        Inputs:
        self :

        Outputs:
        None

        Returns:
        self.comm     : The MPI communicator
        self.mpi_size : The number of the MPI processes in the simulation
        self.mpi_rank : The id of the MPI rank
        self.mpi_stat : A logical variable indicating if we run in parallel or not.
                        If we run on parallel, true; else False

        """
        # Get the MPI information
        self.comm = MPI.COMM_WORLD
        # Get the number of processes
        self.mpi_size = self.comm.Get_size()
        # Get the rank ID
        self.mpi_rank = self.comm.Get_rank()

        # Check if we run serially or in parallel
        if self.mpi_size > 1:
            self.mpi_stat = True
        else:
            self.mpi_stat = False

        # If we run in parallel, gather all the coordinates to be access-
        # ible from all ranks
        if self.mpi_stat:
            # Gather all the coordinates, so as to be accessible from all
            # ranks
            MPITools.gather_coords_across_ranks(self)
            # Calculate the indices where the arrays begin in each rank in
            # relation to the "gathered"/"global" array
            MPITools.get_indices_per_rank(self)

        return


class ForcingTools:
    """

    """

    def update_uv_velocities(self, t, value, item):
        if value == "uv_interpolator":
            ForcingTools.interpolate_uv(
                self,
                t,
                self.uv_x_int,
                self.uv_y_int,
                self.ramp_coeff
            )

        elif value == "uv AMCG interpolator":
            TideTools.set_u_tidal_field(self, t, item)
            TideTools.set_v_tidal_field(self, t, item)


    def configure_uv_interp(self):
        """

        """

        # Load the dictionary relevant to the current forcing
        dict = self.config["uv forcing"]

        if "uv components" in dict:
            # Update the dictionary
            dict = dict["uv components"]
            # If the uv is imported from a file
            if dict["import uv"]:
                ForcingTools.import_uv_from_netCDF(self, dict)

        elif "AMCG forcing" in dict:
            TideTools.configure_amcg_uv_interp(self, dict["AMCG forcing"])

        else:
            mes = ct.br + "Terminating Thetis: The current forcing option " + \
                  "not imported from file has not been programmed yet." + ct.e
            raise SystemExit(mes)
    
    


    def import_uv_from_netCDF(self, dict):
        """

        """
        # If the current forcings are transient
        if dict["transient field"]:
            # Determine the start date of the wind stress
            start_time = TimeTools.calculate_time_since_epoch(
                self, dict["time information"]
            )

            # Get the coordinates to interpolate the uv forcings
            CoordsTools.get_mesh_coords_in_uv_crs(
                self, dict["coordinates information"]
            )

            ForcingTools.import_uv_from_transient_netCDF(self, dict, start_time)



    def import_uv_from_transient_netCDF(self, dict, start_time):
        """

        """
        # Specify the filepath of the wind netcdf
        uv_filepath = dict["nc file"]

        # Determine whether the time information in the netCDF is in the form
        # of seconds since the epoch or in dates (i.e. strings)
        t_type = TimeTools.get_nc_time_info_type(
            self, dict["time information"]
        )

        # Determine whether the wind will be given in as components of a vector
        # or as magnitude and angle
        z1_label, z2_label = IOTools.get_nc_vector_format_and_var_names(
            self, dict["uv information"]
        )

        # If the grid utilised isn't regular structured, i.e. it is rotated or
        # it has not constant spacing or anything else
        if dict["coordinates information"]["rotated grid"]:
            x_data, y_data, t_data, z1_data, z2_data = \
                IOTools.import_transient_irregular_vector_nc(
                    self,
                    uv_filepath,
                    dict["coordinates information"]["x-coords label"],
                    dict["coordinates information"]["y-coords label"],
                    dict["time information"]["label name"],
                    z1_label,
                    z2_label,
                    t_type
                )

        # If the netCDF is structured
        else:
            x_data, y_data, t_data, z1_data, z2_data = \
                IOTools.import_transient_regular_vector_nc(
                    self,
                    uv_filepath,
                    dict["coordinates information"]["x-coords label"],
                    dict["coordinates information"]["y-coords label"],
                    dict["time information"]["label name"],
                    z1_label,
                    z2_label,
                    t_type,
                    dict["coordinates information"]["full coords record"]
                )

        # Process time information
        t_data = TimeTools.process_nc_time_info(
            self,
            t_data,
            start_time,
            dict["time information"]["format"],
            dict["time information"]["timezone"]
        )

        # Convert vector magnitude and direction to xy-components
        if dict["uv information"]["format"] == "magnitude & direction":
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                      "wind vector as magnitude and direction has not been coded" + \
                      " yet." + clrtxt.end
            raise SystemExit(message)

        # Create interpolator
        self.uv_x_int, self.uv_y_int = \
            SetupTools.create_transient_vector_interpolator(
                self,
                x_data, y_data, t_data, z1_data, z2_data,
                dict["coordinates information"]["rotated grid"]
            )

        # Check that we have uv_interpolator in the BC
        ForcingTools.check_vector_uv_update_forcings(self)


        # Interpolate current velocities to mesh grid points
        ForcingTools.interpolate_uv(
            self, 0, self.uv_x_int, self.uv_y_int, self.ramp_coeff
        )


    def check_vector_uv_update_forcings(self):
        """

        """
        # Load dictionary for update forcings
        dict = self.config["update forcings"]

        # Check that we have defined an entry "vectors", since the wind stress
        # is a vector
        if "vectors" not in dict:
            mes = ct.br + "Terminating Thetis: " + ct.e + "If the" + ct.bdc + \
                  " transient field " + ct.e + "in the " + ct.bdc + "uv " + \
                  "forcing " + ct.e + "is " + ct.bdc + "true, " + ct.e + \
                  "then in the entry " + ct.bdc + "update forcings " + ct.e + \
                  "the entry " + ct.bdc + "vectors:{'function_name':'" + \
                  "uv_interpolator'} " + ct.e + "must be included."
            raise SystemExit(mes)

        # Load the dictionary with the vectors in the update forcings key
        dict = dict["vectors"]
        # Get all the items
        items = list(dict.items())

        # Initialise list containing the names of the function to use the uv
        # forcing
        self.uv_func = []

        for item in items:
            if (
                    item[1] == "uv_interpolator" or
                    item[1] == "uv AMCG interpolator"
            ):
                self.uv_func.append(item[0])



        # Check that we have at least one function
        if np.shape(self.uv_func)[0] == 0:
            mes = ct.br + "Terminating Thetis: " + ct.e + "When the field " + \
                  ct.bdc + "transient field " + ct.e + "is " + ct.bdc + \
                  "true" + ct.e + ", then the entry " + ct.bdc + "vectors" + \
                  ct.e + " under the " + ct.bdc + "update forcings" + ct.e + \
                  " needs at least one function with the value " + ct.bdc + \
                  "uv_interpolator" + ct.e + "."
            raise SystemExit(mes)

        return



    def interpolate_uv(self, t, z1_int, z2_int, coeff):
        """

        """
        for func in self.uv_func:
            # Get vector as editable data field
            data = getattr(self, func).dat.data
            data[:,:]= 0.0

            # Create an array for the time dimension
            t = np.ones((np.shape(data[self.bc_ind,:])[0], 1)) * t

            # Concatenate all the coordinates in one array
            points = np.concatenate((t, self.uv_coords), axis=1)

            # Interpolate the x-component
            data[self.bc_ind, 0] = z1_int(points)*coeff
            # Interpolate the y-component
            data[self.bc_ind, 1] = z2_int(points)*coeff

        return



class BCTools(BmiThetis):
    """
    The methods included in this class are:
    * check_functions_existence(self, forc_names)
    * check_function_types(self, user_types)
    * check_user_defined_forcings(self)
    * list_scalar_function(self, t, item)
    * string_scalar_function(self, t, item)
    * update_scalar_forcings(self, t, dict)
    * update_user_defined_forcings(self, t)
    * update_vector_forcings(self, t, dict)
    """
    def __init__(self):
        """

        """
        # The available string options for scalars
        BmiThetis._str_scalar_options = [
            "TPXO forcing", # If the tide forcing is from the TPXO
            "constituents forcing", # If the tide forcing is constant,
            "AMCG forcing", # If the tde consitudents are in AMCG frmat
            "wave_interpolator"  # For WEoC

        ]
        # The available options as function types
        BmiThetis._forcing_func_types = [
            "scalars",
            "vectors"
        ]
        # The available list types of scalar forcings
        BmiThetis._list_scalar_options = [
            "sinusoidal"
        ]

        return


    def check_functions_existence(self, forc_names):
        """
        Check that the functions to be utilised to update the boundary condit-
        ions have been previously defined by the user. If not, terminate Thetis
        Inputs:
        self       :
        forc_names : A 1-D list containing the names of the functions to be us
                     ed to update the forcings

        Outputs:
        None

        Returns:
        None
        """
        # Get the names of the user-defined functions
        func_names = list(self.config["function definition"].keys())

        # Loop through the names of the forcing functions:
        for name in forc_names:
            # If the name is wind_stress_2d
            if name == "wind_stress_2d":
                continue
            if name == "wave":
                continue
            # Check that the forcing function has been defined by the user
            if name not in func_names:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The function under the entry " + \
                    "'update forcings' called '" + clrtxt.bold + \
                    clrtxt.darkcyan + str(name) + clrtxt.end + "' is not " + \
                    "found under the entry 'function definition' in the " + \
                    "Thetis configuration file. The functions, both scalar" + \
                    " and vector ones, under the 'function definition' are" + \
                    clrtxt.bold + clrtxt.darkcyan + str(func_names)[1:-1] + \
                    clrtxt.end
                raise SystemExit(message)

        return


    def check_function_types(self, user_types):
        """
        Check that the types provided under the 'update forcings' entry in the
        json file are matching the approved ones, otherwise terminate Thetis
        Inputs:
        self :
        user_types : A list containing the user-defined function types to be
                     used when updating the forcings

        Outputs:
        None

        Returns:
        None
        """
        # Loop through the functions defined in the json file
        for f_type in user_types:
            # If it isn't one of the available options, terminate Thetis
            if f_type not in self._forcing_func_types:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    "The function type '" + clrtxt.bold + clrtxt.darkcyan + \
                    str(f_type) + clrtxt.end + "' does not match any of " + \
                    "the available function types to be used when updating" + \
                    " the forcings. These are: " + clrtxt.bold + \
                    clrtxt.darkcyan + str(func_types)[1:-1] + clrtxt.end
                raise SystemExit(message)

        return


    def check_user_defined_forcings(self):
        """
        When initialsing Thetis, make sure that the options in update forcings
        are appropriate
        Inputs:
        self :

        Outputs:
        None

        Returns:
        None
        """
        # Initialise BCTools class
        BCTools()

        # Load the dictionary containing the BCs
        dict = self.config["update forcings"]

        # Check the options of the functions to be updated
        BCTools.check_function_types(self, list(dict.keys()))

        ## Check the scalar function forcings
        if "scalars" in dict:
            # Check that the names of the functions for the forcings have been pre-
            # viously defined by the user
            BCTools.check_functions_existence(self, list(dict["scalars"].keys()))

            # Loop through the scalar forcings
            for item in dict["scalars"].items():
                # Check if the value of the key is a list, i.e. it needs more than
                # one number or string to be specified by the user
                if type(item[1]) == list:
                    if item[1][0] not in self._list_scalar_options:
                        mes = ct.br + "Terminating Thetis: " + ct.e + "The " + \
                              "list option " + ct.bdc + str(item[1][0]) + ct.e + \
                              " for " + ct.bdc + "scalar " + ct.e + "forcings " + \
                              "under the entry " + ct.bdc + "update forcings " + \
                              ct.e + "is not included. Please choose one of " + \
                              "the following: " + ct.bdc + \
                              str(self._list_scalar_options)[1:-1] + ct.e + "."
                        raise SystemExit(mes)

                # If the value of the key is a string
                elif type(item[1]) == str:
                    # Check that the string provided is ok, i.e. it corresponds to one of
                    # the forcing options
                    if item[1] not in self._str_scalar_options:
                        mes = ct.br + "Terminating Thetis: " + ct.e + "The " + \
                              "value " + ct.bdc + str(item[1]) + ct.e + " for " + \
                              "scalar forcing under the entry '" + ct.bdc + \
                              "update forcings" + ct.e + "' is not a viable " + \
                              "option. Please choose one of the following: " + \
                              ct.bdc + str(self._str_scalar_options)[1:-1] + \
                              ct.e + "."
                        raise SystemExit(mes)

        ## Check the vector forcings
        if "vectors" in dict:
            # Check that the names of the functions for the forcings have been pre-
            # viously defined by the user
            BCTools.check_functions_existence(self, list(dict["vectors"].keys()))

        return




    def list_scalar_function(self, t, item):
        """
        Assign a scalar function as a forcing boundary condition. The forcing
        parameters and the type of the 'function' to be implemented is given
        through a list
        Inputs:
        self :
        t    : The time to calculate the forcing
        item : The entry (key, value) in the 'update forcings" entry

        Outputs:
        None

        Returns
        """
        # The first item of the list specifies the type of forcing fun-
        # ction, if it is sinusoidal
        if item[1][0] == "sinusoidal":
            # Assign the value to the function
            getattr(self, item[0]).assign(
                item[1][1] * sin(2 * pi / item[1][2] * t)
            )
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                clrtxt.end + "A scalar forcing function under the entry " + \
                "'update forcings' containing a list as its entry should " + \
                "have as its first item '" + clrtxt.bold + clrtxt.darkcyan + \
                "sinusoidal" + clrtxt.end + "' and not '" + clrtxt.bold + \
                clrtxt.darkcyan + str(itwm[1][0]) + clrtxt.end + "'."
            raise SystemExit(message)

        return


    def string_scalar_function(self, t, item):
        """

        """
        # If the value is "TPXO"
        if item[1] == "TPXO forcing" or item[1] == "AMCG forcing":
            TideTools.set_tidal_field(self, t, item)

        elif item[1] == "constituents forcing":
            TideTools.update_tide_time_function(self, t, item)

        elif item[1] == "wave_interpolator":
            WEoCTools.interpolate_wave_char(self, t)

        return



    def update_scalar_forcings(self, t, dict):
        """
        Update the scalar functions that represent the forcings for the Thetis
        simulation, based on the type of the value, list, string, number
        Inputs:
        self :
        t    : The time to calculate the forcing
        dict : A dictionary containing the scalar function forcing parameters

        Outputs:
        None

        Returns:
        None
        """
        # Loop through the forcings
        for item in dict.items():
            # Check if the value of the key is a list, i.e. it needs more than
            # one number or string to be specified by the user
            if type(item[1]) == list:
                BCTools.list_scalar_function(self, t, item)

            # If the value of the key is a string
            elif type(item[1]) == str:
                BCTools.string_scalar_function(self, t, item)

            # If the value of the key is a number
            elif ( type(item[1]) == float or type(item[1]) == int ):
                # Assign the constant value
                getattr(self, item[0]).assign(item[1])

            # If the value of the key is anything else, terminate Thetis
            else:
                message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                    clrtxt.end + "The value of the key '" + clrtxt.bold + \
                    clrtxt.darkcyan + str(item[0]) + clrtxt.end + "' is " + \
                    "not a " + clrtxt.bold + clrtxt.darkcyan + "list, " + \
                    "string, " + clrtxt.end + "or a " + clrtxt.bold + \
                    clrtxt.darkcyan + "number" + clrtxt.end + "."
                raise SystemExit(message)

        return


    def update_user_defined_forcings(self, t):
        """

        """
        # Load the dictionary containing the BCs
        dict = self.config["update forcings"]

        # If the key "scalars" exists, update the scalar functions
        if "scalars" in dict:
            BCTools.update_scalar_forcings(self, t, dict["scalars"])

        # If the key "vectors" exists, update the vector functions
        if "vectors" in dict:
            BCTools.update_vector_forcings(self, t, dict["vectors"])

        return


    def update_vector_forcings(self, t, dict):
        """
        Update the vector functions that represent the forcings for the Thetis
        simulation. The only option is through a list that is passed as a
        vector to the vector function
        Inputs:
        self :
        t    : The time to calculate the forcing
        dict : A dictionary containing the vector forcing parameters

        Outputs:
        None

        Returns:
        None
        """
        # Loop through the forcings
        for item in dict.items():
            # If the vector is wind stress
            if item[0] == "wind_stress_2d":
                WindTools.update_wind_stress(self, t, dict["wind_stress_2d"])

            # Else if the function is one of the velocity forcings
            elif item[0] in self.uv_func:
                ForcingTools.update_uv_velocities(self, t, dict[item[0]], item[0])

            # For the other vectors
            else:
                getattr(self, item[0]).assign(as_vector(item[1]))

        return



class TimeTools(BmiThetis):
    """
    This class contains the following methods:
    * calculate_time_since_epoch(self, dict)
    * create_UTC_date(self, dict)
    * get_datetime(self, dict)
    * get_nc_time_info_type(self, dict)
    * process_nc_time_info(self, t_data, start_time, t_format)
    """

    def calculate_time_since_epoch(self, dict):
        """
        Calculate the time since the epoch time, i.e. 01/01/1970, for the time
        described in the provided dictionary
        Inputs:
        self :
        dict : A dictionary containing the parameters for TPXO forcing

        Outputs:
        None

        Returns:
        time : The time since epoch of the specified date [s]
        """
        # Create datetime from dictionary information
        date = TimeTools.get_datetime(self, dict)

        # Conevrt datetime to time since epoch
        time = date.timestamp()

        return time


    def create_UTC_date(self, dict):
        """
        Create a datetime variable in the UTC timezone
        Inputs:
        self :
        dict : A dictionary containing the parameters for TPXO forcing

        Outputs:
        None

        Returns:
        date : The datetime of the provided information
        """
        import datetime

        # Calculate the start datetime
        date = datetime.datetime(
            dict["start year"],
            dict["start month"],
            dict["start day"],
            dict["start hour"],
            dict["start minute"],
            tzinfo = datetime.timezone.utc
        )

        return date


    def get_datetime(self, dict):
        """
        Get a datetime 'object' depending on the timezone
        Inputs:
        self :
        dict : A dictionary contaning the datetime information

        Outputs:
        None

        Returns:
        date : A datetime
        """
        # If the timezone is UTC:
        if dict["timezone"] == "UTC":
            date = TimeTools.create_UTC_date(self, dict)

        # if it isn't
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                "The creation of a datetime with a timezone other than " + \
                "UTC has not been coded yet"

        return date


    def get_nc_time_info_type(self, dict):
        """
        Get the type of the time information available in the netCDF.
        Inputs:
        self:
        dict : Dictionary containing the time information of the netCDF

        Outputs:
        None

        Returns:
        t_type : A string containing the type to be loaded from the netcdf
        """
        # Get the format of the time information
        t_format = dict["format"]

        # If the time is given in seconds since the epoch time
        if t_format == "Time since epoch [s]":
            t_type = 'float64'

        # If the time is given in datetime type as DD/MM/YYYY;hh:mm:ss
        elif t_format == 'Datetime DD/MM/YYYY;hh:mm:ss':
            t_type = 'str'

        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                "No other time information type has been coded yet apart " + \
                "from " + clrtxt.bold + clrtxt.darkcyan + "'Time since " + \
                "epoch [s]'" + clrtxt.end + " and " + clrtxt.bold + \
                clrtxt.darkcyan + "'Datetime DD/MM/YYYY;hh:mm:ss'" + \
                clrtxt.end
            raise SystemExit(message)

        return t_type


    def process_nc_time_info(self, t_data, start_time, t_format, t_zone):
        """
        Process the time information loaded from a netCDF file. This means con-
        vert it to time since epoch if it is in another format and then calcul-
        ate the relative times based on the start time of simulation/forcing.
        Inputs:
        self       :
        t_data     : 1-D array containing the time information from the netCDF
        start_time : The start time of the simulation since epoch [s]
        t_format   : The format of the time information of t_data
        t_zone     : The timezone utilised in the time information

        Outputs:
        None

        Returns:
        t_data : 1-D array containing the relative times of the netcdf with re-
                 gards the start_time [s]
        """

        # Convert datetimes to time since epoch
        if t_format == 'Datetime DD/MM/YYYY;hh:mm:ss':
            message = clrtxt.bold + clrtxt.red + "Termintaing Thetis: " + \
                "The procession of the nc time information with datetime " + \
                "format is not coded yet" + clrtxt.end
            raise SystemExit

        # Convert the time since epoch to timezone
        t_data = t_data - start_time

        return t_data



class WindTools(BmiThetis):
    """
    This class contains the following methods:
    * __init__(self)
    * check_vector_wind_update_forcings(self)
    * check_wind_drag_formulation(self, dict)
    * import_wind_from_netCDF(self, dict)
    * import_wind_from_transient_netCDF(self, dict, start_time)
    * interpolate_wind_stress( self, t, z1_int, z2_int)
    * interpolate_wind_velocities(self, t, z1_int, z2_int)
    * update_wind_stress(self, t)
    * wind_field(self)

    """

    def __init__(self):
        """
        Create a BMI-refactored 2-D Thetis model that is ready for initialisa-
        tion
        """
        # The available methods/formulations to calculate the wind drag coeffi-
        # cient C_D
        BmiThetis._wind_methods = [
            "LargePond1981",
            "SmithBanke1975",
            "Wu1969",
            "LargeYeager2009"
        ]

        # The available options for the update forcings for the wind_stress_2d
        BmiThetis._wind_vector_update_forc_val = [
            "interpolator"
        ]
        return

    def check_vector_wind_update_forcings(self):
        """

        """
        # Load dictionary for update forcings
        dict = self.config["update forcings"]

        # Check that we have defined an entry "vectors", since the wind stress
        # is a vector
        if "vectors" not in dict:
            mes = ct.br + "Terminating Thetis: " + ct.e + "If the" + ct.bdc + \
                  " transient field " + ct.e + "in the " + ct.bdc + "wind " + \
                  "forcing " + ct.e + "is " + ct.bdc + "true, " + ct.e + \
                  "then in the entry " + ct.bdc + "update forcings " + ct.e + \
                  "the entry " + ct.bdc + "vectors:{'wind_stress_2d':'" + \
                  "interpolator'} " + ct.e + "must be included."
            raise SystemExit(mes)

        # Check that the wind_stress_2d is included in the vectors dictionary
        if "wind_stress_2d" not in dict["vectors"]:
            mes = ct.br + "Terminating Thetis: " + ct.e + "The item " + \
                  ct.bdc + "'wind_stress_2d':'interpolator' " + ct.e + \
                  "must be included under the entry " + ct.bdc + "vectors " + \
                  ct.e + "located under the " + ct.bdc + "update forcings " + \
                  ct.e + "to update the wind stress"
            raise SystemExit(mes)

        # Check that the wind_stress_2d key has the value "interpolator"
        value = dict["vectors"]["wind_stress_2d"]
        avail_values = self._wind_vector_update_forc_val
        if value not in avail_values:
            mes = ct.br + "Terminating Thetis: " + ct.e + "The value of " + \
                  "the entry " + ct.bdc + "wind_stress_2d " + ct.e + \
                  "under the entry " + ct.bdc + "vectors " + ct.e + "of " + \
                  "the entry " + ct.bdc + "update forcings " + ct.e + \
                  "must be one of the following: " + ct.bdc + \
                  str(avail_values)[1:-1] + ct.e + "."
            raise SystemExit(mes)





    def check_wind_drag_formulation(self, dict):
        """
        Check that the user has specified a wind drag formulation to be used
        and that the one specified is one of the available ones described in
        self._wind_methods
        Inputs:
        self :
        dict : A dictionary containing the wind information

        Outputs:
        None

        Returns:
        None
        """
        # Check that the key 'formulation' exists in the dictionary
        if "formulation" not in dict:
            mes = ct.br + "Terminating Thetis: " + ct.e + "If the entry '" + \
                  "field" + ct.e + "' under the entry '" + ct.bdc + "wind " + \
                  "information" + ct.e + "' is '" + ct.bdc + "wind " + \
                  "velocities" + ct.e + "', then the entry '" + ct.bdc + \
                  "formulation" + ct.e + "' should also exist under the " + \
                  "latter. The available formulations are: " + ct.bdc + \
                   str(self._wind_methods)[1:-1] + ct.e
            raise SystemExit(mes)

        if dict["formulation"] not in self._wind_methods:
            mes = ct.br + "Terminating Thetis:  " + ct.e + "The wind drag " + \
                  "formulation " + ct.bdc + str(dict["formulation"]) + ct.e + \
                  " is not available. Please choose one of the following: " + \
                  ct.bdc + str(self._wind_methods)[1:-1] + ct.e
            raise SystemExit(mes)











    def import_wind_from_netCDF(self, dict):
        """

        """
        # If the wind forcings are transient with time
        if dict["transient field"]:
            # Determine the start date of the wind stress
            start_time = TimeTools.calculate_time_since_epoch(
                self, dict["time information"]
            )

            # Get the coordinates to interpolate the wind forcings
            CoordsTools.get_mesh_coords_in_wind_crs(
                self, dict["coordinates information"]
            )

            WindTools.import_wind_from_transient_netCDF(self, dict, start_time)

            # Check that we have included the wind stresses in the update forc-
            #ings
            WindTools.check_vector_wind_update_forcings(self)

        # If the wind forcings are constant with time
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                "A stationary wind forcing imported from netcdf has not " + \
                "coded yet"
            raise SystemExit(message)


    def import_wind_from_transient_netCDF(self, dict, start_time):
        """

        """
        # Specify the filepath of the wind netcdf
        wind_filepath = dict["nc file"]

        # Determine whether the time information in the netCDF is in the form
        # of seconds since the epoch or in dates (i.e. strings)
        t_type = TimeTools.get_nc_time_info_type(
            self, dict["time information"]
        )

        # Determine whether the wind will be given in as components of a vector
        # or as magnitude and angle
        z1_label, z2_label = IOTools.get_nc_vector_format_and_var_names(
            self, dict["wind information"]
        )

        # If the grid utilised isn't regular structured, i.e. it is rotated or
        # it has not constant spacing or anything else
        if dict["coordinates information"]["rotated grid"]:
            x_data, y_data, t_data, z1_data, z2_data = \
                IOTools.import_transient_irregular_vector_nc(
                    self,
                    wind_filepath,
                    dict["coordinates information"]["x-coords label"],
                    dict["coordinates information"]["y-coords label"],
                    dict["time information"]["label name"],
                    z1_label,
                    z2_label,
                    t_type
            )

        # If the netCDF is structured
        else:
            x_data, y_data, t_data, z1_data, z2_data = \
                IOTools.import_transient_regular_vector_nc(
                    self,
                    wind_filepath,
                    dict["coordinates information"]["x-coords label"],
                    dict["coordinates information"]["y-coords label"],
                    dict["time information"]["label name"],
                    z1_label,
                    z2_label,
                    t_type,
                    dict["coordinates information"]["full coords record"]
            )

        # Process time information
        t_data = TimeTools.process_nc_time_info(
            self,
            t_data,
            start_time,
            dict["time information"]["format"],
            dict["time information"]["timezone"]
        )

        # Convert vector magnitude and direction to xy-components
        if dict["wind information"]["format"] == "magnitude & direction":
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                "wind vector as magnitude and direction has not been coded" + \
                " yet." + clrtxt.end
            raise SystemExit(message)

        # Create interpolator
        self.wind_x_int, self.wind_y_int = \
            SetupTools.create_transient_vector_interpolator(
                self,
                x_data, y_data, t_data, z1_data, z2_data,
                dict["coordinates information"]["rotated grid"]
        )

        # Define Vector Function for wind velocities
        self.wind_uv = Function(self.P1v_2d, name="wind_uv_2d")
        # Define Function for wind velocity magnitude
        self.wind_uv_magn = Function(self.P1_2d, name="wind_uv_magn")

        # Interpolate wind stress to mesh grid points
        WindTools.interpolate_wind_stress(
            self, 0, self.wind_x_int, self.wind_y_int, dict["wind information"]
        )
        # if dict["wind information"]["field"] == "wind velocities":



    def interpolate_wind_stress(self, t, z1_int, z2_int, dict):
        """

        """
        from thetis.forcing import compute_wind_stress
        if dict["field"] == "wind velocities":


            # Vector components as editable fields
            ws_x = self.wind_stress_2d.dat.data[:,0]
            ws_y = self.wind_stress_2d.dat.data[:,1]

            # Interpolate wind velocities to mesh point
            WindTools.interpolate_wind_velocities(self, t, z1_int, z2_int)
            # Check that the fomulation for the calculation of the wind drag
            # coefficient is valid
            WindTools.check_wind_drag_formulation(self, dict)
            # WindTools.calculate_wind_stress(self, dict)
            # Calculate wind stress
            ws_x, ws_y = compute_wind_stress(
                self.wind_uv.dat.data[:,0],
                self.wind_uv.dat.data[:,1],
                method = dict["formulation"]
            )

            self.wind_stress_2d.dat.data[:, 0] = ws_x
            self.wind_stress_2d.dat.data[:, 1] = ws_y


        # # Vector as editable field
        # data = self.wind_stress_2d .dat.data
        # print(data[:,0])
        # # Create an array for the time dimension
        # t = np.ones( (np.shape(data)[0],1)) * t
        #
        # # Concatenate all the coordinates in one array
        # points = np.concatenate( (t, self.wind_coords), axis=1)
        #
        # # Interpolate the x-component
        # data[:,0] = z1_int(points)
        # # Interpolate the y-component
        # data[:,1] = z2_int(points)
        #
        # # Check if the x,y-components are wind stress or wind velocities
        # if dict["field"] == "wind velocities":
        #     WindTools.calculate_wind_stress(self, dict)



    def interpolate_wind_velocities(self, t, z1_int, z2_int):
        """
        Interpolate the wind velocities at time t across the mesh grid points,
        utilising the provided interpolators for each component
        Inputs:
        self :
        t    :
        z1_int :
        z2_int :

        Outputs:
        None

        Returns:
        self.wind_uv : The wind velocities at the mesh grid points
        """
        # Vector as editable field
        data = self.wind_uv.dat.data

        # Create an array for the time dimension
        t = np.ones((np.shape(data)[0], 1)) * t

        # Concatenate all the coordinates in one array
        points = np.concatenate((t, self.wind_coords), axis=1)

        # Interpolate the x-component
        data[:, 0] = z1_int(points)
        # Interpolate the y-component
        data[:, 1] = z2_int(points)


        return


    def update_wind_stress(self, t, value):
        """
        Update the wind stress field depending on the value of the key
        "wind_stress_2d" under the entry "vectors" in the "update forcings".
        Inputs:
        self  :
        t     : The time of the simulation
        value : The value of the wins_stress_2d key in the update_forcings

        Outputs:
        None

        Returns:
        None
        """
        if value=="interpolator":
            WindTools.interpolate_wind_stress(
                self,
                t,
                self.wind_x_int,
                self.wind_y_int,
                self.config["wind forcing"]["wind information"]
            )

        return


    def wind_field(self):
        """

        """
        from netCDF4 import Dataset
        # Initialise WInd Class?
        WindTools()

        # Define VectorFunction
        self.wind_stress_2d = Function(self.P1v_2d, name="wind_stress")

        # Load the dictionary relevant to the wind forcing
        dict = self.config["wind forcing"]

        # If the wind is imported from a file
        if dict["import wind"]:
            WindTools.import_wind_from_netCDF(self, dict)

        # If it isn't, it should be constant
        else:
            message = clrtxt.bold + clrtxt.red + "Terminating Thetis: " + \
                "Wind forcing not imported from a netCDF has not been " + \
                "coded yet" + clrtxt.end
            raise SystemExit(message)



class FinaliseTools(BmiThetis):
    """
    This class contains the following methods:
    * __init__(self)
    * check_finalise_params(self)
    * check_time_info(self)
    """
    def __init__(self):
        """

        """
        BmiThetis._final_options = [
            "output directory",
            "time information"
        ]

        return


    def check_finalise_params(self):
        """
        Check that the user has included a dictionary containing the parameters
        for when Thetis is exiting the time-loop. Also confirm that the param-
        eters specified are appropriate. Otherwise, terminate the model
        Inputs:
        self :

        Outputs:
        None

        Returns:
        None
        """
        # Initialise FinaliseTools class
        FinaliseTools()

        # Check that the user has defined the entry "finalize parameters"
        if "finalize parameters" not in self.config:
            mes = ct.br + "Terminating Thetis: " + ct.e + "The entry " + \
                  ct.bdc + "finalize parameters " + ct.e + "must exist in " + \
                  "the configuration file. Please add it."
            raise SystemExit(mes)

        # Load the dictionary containing the exiting parameters
        dict = self.config["finalize parameters"]

        # Loop through the kesy of the dictionary
        for key in list(dict.keys()):
            # Check if the keys are the appropriate ones
            if key not in self._final_options:
                mes = ct.br + "Terminating Thetis: " + ct.e + "They key " + \
                      ct.bdc + str(key) + ct.e + " is not a viable entry. " + \
                      "Please choose one of the following: " + ct.bdc + \
                    str(self._final_options)[1:-1] + ct.e + "."
                raise SystemExit(mes)

        return


    def check_time_info(self):
        """

        """
        # Load dictionary containing the parameters for finalisation
        dict = self.config["finalize parameters"]

        # Check if the "time information" exists
        if "time information" in dict:
            time_info = dict['time information']

        # If it doesn't exist
        else:
            time_info = None

        return time_info


class WEoCTools(BmiThetis):
    """

    """
    def __init__(self):
        BmiThetis._format_options = [
            "characteristics"
        ]
        BmiThetis._wave_labels_options = [
            "hwave label",
            "wlen label",
            "dir label",
            "qb label"
        ]

    def configure_weoc_forcing(self):
        """
        Setup an interpolator (most likely) to account for the Wave Effects on
        Currents.

        Returns:

        """
        # Initialise WEoCTools class
        WEoCTools()

        # Load dictionary rearding the Wave Effects on Currents:
        dict = self.config["WEoC"]

        # Check that the format is approprita
        if dict["format"] not in self._format_options:
            mes = f"{ct.br}Terminating Thetis: {ct.e}The option {ct.bdc}" + \
                  f"{dict['format']}{ct.e} is not a viable options for " + \
                  f"WEoC forcing. Please choose one fo the following: " + \
                  f"{ct.bdc}{str(self._format_options)[1:-1]}{ct.e}"
            raise SystemExit(mes)

        if dict["format"] == "characteristics":
            WEoCTools.create_wave_char_interp(self, dict)


    def create_wave_char_interp(self, dict):
        """

        Args:
            dict:

        Returns:

        """
        # Setup the functions for H, Dir, WLen and Qb
        self.hwave = Function(self.P1_2d, name="wave_height")
        self.dir = Function(self.P1_2d, name="mean_wave_direction")
        self.wlen = Function(self.P1_2d, name="mean_wavelength")
        self.qb = Function(self.P1_2d, name="wave_breaking_percentage")


        # Check that we have all the labels
        for label in BmiThetis._wave_labels_options:
            if label not in dict:
                mes = f"{ct.br}Terminating Thetis:{ct.e} The entry " + \
                      f"{ct.bdc}{label}{ct.e} is not in the entry " + \
                      f"{ct.bdc}WEoC{ct.e}. Please include it"
                raise SystemExit(mes)

        # If the Wave Effects on Currents are transient
        if dict["transient field"]:
            # Determine the start date of the wind stress
            start_time = TimeTools.calculate_time_since_epoch(
                self, dict["time information"]
            )

            # Get the coordinates to interpolate the uv forcings
            CoordsTools.get_mesh_coords_in_wave_char_crs(
                self, dict["coordinates information"]
            )

            WEoCTools.import_wave_char_from_transient_netCDF(
                self, dict, start_time
            )



    def import_wave_char_from_transient_netCDF(self, dict, start_time):
        """
        Import the specified parameter form the netcdf
        Args:
            dict:
            start_time:
            label:

        Returns:

        """
        # Specify the file to be opened
        wave_fpath = dict["nc file"]

        # Determine whether the time information in the netCDF is in the form
        # of seconds since the epoch or in dates (i.e. strings)
        t_type = TimeTools.get_nc_time_info_type(
            self, dict["time information"]
        )

        for label in self._wave_labels_options:
            # If the grid utilised isn't regular structured, i.e. it is rotated or
            # it has not constant spacing or anything else
            if dict["coordinates information"]["rotated grid"]:
                x_data, y_data, t_data, z_data = \
                    IOTools.import_transient_irregular_scalar_nc(
                        self,
                        wave_fpath,
                        dict["coordinates information"]["x-coords label"],
                        dict["coordinates information"]["y-coords label"],
                        dict["time information"]["label name"],
                        dict[label],
                        t_type
                    )

            # If the netCDF is structured
            else:
                x_data, y_data, t_data, z_data = \
                    IOTools.import_transient_regular_scalar_nc(
                        self,
                        wave_fpath,
                        dict["coordinates information"]["x-coords label"],
                        dict["coordinates information"]["y-coords label"],
                        dict["time information"]["label name"],
                        dict[label],
                        t_type,
                        dict["coordinates information"]["full coords record"]
                    )

            # Process time information
            t_data = TimeTools.process_nc_time_info(
                self,
                t_data,
                start_time,
                dict["time information"]["format"],
                dict["time information"]["timezone"]
            )

            if label == "hwave label":
                # Create hwave interpolator
                self.hwave_int = \
                    SetupTools.create_transient_scalar_interpolator(
                        self,
                        x_data, y_data, t_data, z_data,
                        dict["coordinates information"]["rotated grid"]
                    )
            elif label == "wlen label":
                self.wlen_int = \
                    SetupTools.create_transient_scalar_interpolator(
                        self,
                        x_data, y_data, t_data, z_data,
                        dict["coordinates information"]["rotated grid"]
                    )
            elif label == "dir label":
                self.dir_int = \
                    SetupTools.create_transient_scalar_interpolator(
                        self,
                        x_data, y_data, t_data, z_data,
                        dict["coordinates information"]["rotated grid"]
                    )

            elif label == "qb label":
                self.qb_int = \
                    SetupTools.create_transient_scalar_interpolator(
                        self,
                        x_data, y_data, t_data, z_data,
                        dict["coordinates information"]["rotated grid"]
                    )

        # Check that we have the necessary interpolators
        # TO-DO

        # Interpolate the wave characteristics
        WEoCTools.interpolate_wave_char(self, 0)

        return



    def interpolate_wave_char(self, t ):
        """

        Args:
            t:
            z_int:
            coeff:

        Returns:

        """
        # Create an array for the time dimension
        t = np.ones((np.shape(self.hwave.dat.data)[0], 1)) * t

        # Concatenate all the coordinates in one array
        points = np.concatenate((t, self.wave_coords), axis=1)


        # Interpolate
        self.hwave.dat.data[:] = self.hwave_int(points)
        self.dir.dat.data[:] = self.dir_int(points)
        self.wlen.dat.data[:] = self.wlen_int(points)
        self.qb.dat.data[:] = self.qb_int(points)

        return
















