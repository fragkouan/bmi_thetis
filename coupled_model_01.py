"""
Script to run the coupled model
"""
# bmithetis should always be imported before bmi_SWAN, else we have issues
from bmithetis import BmiThetis
from bmi_swan import Swan
from coupling_checks import *
from datetime import timedelta
import datetime
from math import ceil
import numpy as np
from thetis import *
import time

# Start time
t1 = time.time()
# Print datetime

# SWAN configuration file
config_file = "./westray_01.swn"

# Initialize Thetis
thetis = BmiThetis()
thetis.initialize('./inputs/thetis_01_configuration.json')



# Initialize SWAN
swan = Swan()
swan.initialize(config_file)

# Check the coupled status for both models
coupled_stat = check_coupling_status(
    thetis.coupl_stat,
    swan.get_value(
        'coupling_status',
        np.empty(1,swan.get_var_type('coupling_status'))
    )[0]
)


if coupled_stat!='No coupling':
    # Check the coupling interval
    check_coupling_dt(
        thetis.coupl_dt,
        swan.get_value(
            'coupling_timestep',
            np.empty(1, swan.get_var_type('coupling_timestep'))
        )[0]
    )
    # Check the duration of the simulation
    check_simulation_duration(
        thetis.t_end-thetis.t_start,
        swan.get_value(
            "simulation_duration",
            np.empty(1, swan.get_var_type('simulation_duration'))
        )[0]
    )
else:
    raise SystemExit("No coupling was selected")

# Get the necessary parameters, liek start time, end time and coupling timestep
# to create ana rray of the appropriate times
times = np.arange(
    thetis.t_start ,
    thetis.t_end + thetis.coupl_dt,
    thetis.coupl_dt
)

if ((coupled_stat=='Fully coupled') or (coupled_stat=='SWAN-to-Thetis')):
    # Create the arrays for H, Dir, Wlen, Qb
    hwave_data = np.empty(
        swan.get_grid_size(
            swan.get_var_grid("sea_surface_water_wave__height")
        ),
        dtype=swan.get_var_type("sea_surface_water_wave__height")
    )
    dir_data = np.empty(
        swan.get_grid_size(
            swan.get_var_grid("sea_surface_water_wave__direction")
        ),
        dtype=swan.get_var_type("sea_surface_water_wave__direction")
    )
    wlen_data = np.empty(
        swan.get_grid_size(
            swan.get_var_grid("sea_surface_water_wave__wavelength")
        ),
        dtype=swan.get_var_type("sea_surface_water_wave__wavelength")
    )
    qb_data = np.empty(
        swan.get_grid_size(
            swan.get_var_grid("sea_surface_water_wave__breaking_fraction")
        ),
        dtype=swan.get_var_type("sea_surface_water_wave__breaking_fraction")
    )

if ((coupled_stat=='Fully coupled') or (coupled_stat=='Thetis-to-SWAN')):
    # Create the arrays for elev, u, and v
    wlevl_data = np.empty(
        swan.get_grid_size(
            swan.get_var_grid("sea_water_surface__elevation")
        )
    )
    uxb_data = np.empty(
        swan.get_grid_size(
            swan.get_var_grid("sea_water_flow__x_component_of_velocity")
        )
    )
    uyb_data = np.empty(
        swan.get_grid_size(
            swan.get_var_grid("sea_water_flow__y_component_of_velocity")
        )
    )
# wlevl_data = thetis.get_value('sea_water_surface__elevation')
# uxb_data = thetis.get_value('sea_water_flow__x_component_of_velocity')
# uyb_data = thetis.get_value('sea_water_flow__y_component_of_velocity')

# Loop through times, i.e. Run the simulation
for t in times:
    print_output("Time = " + str(t))
    if t==0:
        if coupled_stat=='Fully coupled' or coupled_stat=='Thetis-to-SWAN':
            wlevl_data = thetis.get_value('sea_water_surface__elevation')
            uxb_data = thetis.get_value('sea_water_flow__x_component_of_velocity')
            uyb_data = thetis.get_value('sea_water_flow__y_component_of_velocity')
            swan.set_value(
                'sea_water_surface__elevation',
                wlevl_data.astype(
                    swan.get_var_type(
                        "sea_water_surface__elevation"
                    )
                )
            )
            swan.set_value(
                'sea_water_flow__x_component_of_velocity',
                uxb_data.astype(
                    swan.get_var_type(
                        "sea_water_flow__x_component_of_velocity"
                    )
                )
            )
            swan.set_value(
                'sea_water_flow__y_component_of_velocity',
                uyb_data.astype(
                    swan.get_var_type(
                        "sea_water_flow__x_component_of_velocity"
                    )
                )
            )
        swan.update_until(t)
        # Get the necessary infromation from SWAN
        if coupled_stat=='Fully coupled' or coupled_stat=="SWAN-to-Thetis":
            swan.get_value("sea_surface_water_wave__height", hwave_data)
            swan.get_value("sea_surface_water_wave__direction", dir_data)
            swan.get_value("sea_surface_water_wave__wavelength", wlen_data)
            swan.get_value("sea_surface_water_wave__breaking_fraction", qb_data)

    else:
        if coupled_stat=='Fully coupled' or coupled_stat=='Thetis-to-SWAN':
            wlevl_data = thetis.get_value('sea_water_surface__elevation')
            uxb_data = thetis.get_value('sea_water_flow__x_component_of_velocity')
            uyb_data = thetis.get_value('sea_water_flow__y_component_of_velocity')
            swan.set_value(
                'sea_water_surface__elevation',
                wlevl_data.astype(
                    swan.get_var_type(
                        "sea_water_surface__elevation"
                    )
                )
            )
            swan.set_value(
                'sea_water_flow__x_component_of_velocity',
                uxb_data.astype(
                    swan.get_var_type(
                        "sea_water_flow__x_component_of_velocity"
                    )
                )
            )
            swan.set_value(
                'sea_water_flow__y_component_of_velocity',
                uyb_data.astype(
                    swan.get_var_type(
                        "sea_water_flow__x_component_of_velocity"
                    )
                )
            )
        swan.update_until(t)
        # Get the necessary infromation from SWAN
        if coupled_stat=='Fully coupled' or coupled_stat=="SWAN-to-Thetis":
            swan.get_value("sea_surface_water_wave__height", hwave_data)
            swan.get_value("sea_surface_water_wave__direction", dir_data)
            swan.get_value("sea_surface_water_wave__wavelength", wlen_data)
            swan.get_value("sea_surface_water_wave__breaking_fraction", qb_data)

        # Pass SWAN fields to Thetis
        if coupled_stat=='Fully coupled' or coupled_stat=="SWAN-to-Thetis":
            # if t > thetis.coupl_dt:
            thetis.set_value("sea_surface_water_wave__height", hwave_data)
            thetis.set_value("sea_surface_water_wave__direction", dir_data)
            thetis.set_value("sea_surface_water_wave__wavelength", wlen_data)
            thetis.set_value("sea_surface_water_wave__breaking_fraction", qb_data)

        thetis.update_until(t)

        # Get the necessary fields' values from Thetis
        if coupled_stat=='Fully coupled' or coupled_stat=='Thetis-to-SWAN':
            wlevl_data = thetis.get_value('sea_water_surface__elevation')
            # print(f"eta bmi:[{np.amin(wlevl_data):.3f}, {np.amax(wlevl_data):.3f}]")
            uxb_data = thetis.get_value('sea_water_flow__x_component_of_velocity')
            uyb_data = thetis.get_value('sea_water_flow__y_component_of_velocity')
            # print(f"eta=[{min(wlevl_data):.3f},{max(wlevl_data):.3f}] m")
            swan.set_value(
                'sea_water_surface__elevation',
                wlevl_data.astype(
                    swan.get_var_type(
                        "sea_water_surface__elevation"
                    )
                )
            )
            swan.set_value(
                'sea_water_flow__x_component_of_velocity',
                uxb_data.astype(
                    swan.get_var_type(
                        "sea_water_flow__x_component_of_velocity"
                    )
                )
            )
            swan.set_value(
                'sea_water_flow__y_component_of_velocity',
                uyb_data.astype(
                    swan.get_var_type(
                        "sea_water_flow__x_component_of_velocity"
                    )
                )
            )



# Finalise model
thetis.finalize()
# Tidy-up and "close" SWAN
swan.finalize()


# End time
t2 = time.time()

# Print time
td_str =str(timedelta(seconds=t2-t1))
x = td_str.split(':')
x[2] = ceil(float(x[2]))
print("Simulation duration : ",x[0],"Hours",x[1],"Minutes",x[2],"Seconds")



