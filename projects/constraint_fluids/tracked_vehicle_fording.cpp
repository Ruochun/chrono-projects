// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
//
// =============================================================================

#include <iostream>

// Chrono::Engine header files
#include "chrono/ChConfig.h"
#include "chrono/core/ChFileutils.h"
#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsInputOutput.h"

// Chrono::Parallel header files
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"
#include "chrono_parallel/collision/ChNarrowphaseRUtils.h"

// Chrono::Parallel OpenGL header files
//#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

// Chrono utility header files
#include "utils/ChUtilsGeometry.h"
#include "utils/ChUtilsCreators.h"
#include "utils/ChUtilsGenerators.h"
#include "utils/ChUtilsInputOutput.h"

// Chrono vehicle header files
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/utils/ChSpeedController.h"
// M113 model header files
#include "chrono_models/vehicle/m113/M113_SimplePowertrain.h"
#include "chrono_models/vehicle/m113/M113_Vehicle.h"

#include "fording_setup.h"


using namespace chrono;
using namespace chrono::collision;
using namespace chrono::vehicle;
using namespace chrono::vehicle::m113;

using std::cout;
using std::endl;

// =============================================================================
// USER SETTINGS
// =============================================================================


// -----------------------------------------------------------------------------
// Specification of the terrain
// -----------------------------------------------------------------------------

enum TerrainType { RIGID_TERRAIN, GRANULAR_TERRAIN };

// Type of terrain
TerrainType terrain_type = RIGID_TERRAIN;

// Control visibility of containing bin walls
bool visible_walls = false;

// Dimensions
double hdimX = 5.5; //// 2.5;
double hdimY = 2.5;
double hdimZ = 0.5;
double hthick = 0.25;

// Parameters for granular material
int Id_g = 100;
double r_g = 0.02;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

float mu_g = 0.8f;

unsigned int num_particles = 100; //// 40000;

// -----------------------------------------------------------------------------
// Specification of the vehicle model
// -----------------------------------------------------------------------------
// Initial vehicle position and orientation

// Simple powertrain model
std::string simplepowertrain_file("generic/powertrain/SimplePowertrain.json");

ChSpeedController m_speedPIDVehicle;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Total simulation duration.
double time_end = 20;

// Duration of the "hold time" (vehicle chassis fixed and no driver inputs).
// This can be used to allow the granular material to settle.
double time_hold = 0.1;

// Solver parameters
double time_step = 1e-3;  // 2e-4;

double tolerance = 0.001;

int max_iteration_bilateral = 1000;  // 1000;
int max_iteration_normal = 0;
int max_iteration_sliding = 50;  // 2000;
int max_iteration_spinning = 0;

float contact_recovery_speed = -1;

// Periodically monitor maximum bilateral constraint violation
bool monitor_bilaterals = false;
int bilateral_frame_interval = 100;

int out_fps = 60;

double target_speed = 2;
// =============================================================================

double CreateParticles(ChSystem* system) {
	// Create a material
	auto mat_g = std::make_shared<ChMaterialSurface>();
	mat_g->SetFriction(mu_g);

	// Create a particle generator and a mixture entirely made out of spheres
	utils::Generator gen(system);
	std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
	m1->setDefaultMaterial(mat_g);
	m1->setDefaultDensity(rho_g);
	m1->setDefaultSize(r_g);

	// Set starting value for body identifiers
	gen.setBodyIdentifier(Id_g);

	// Create particles in layers until reaching the desired number of particles
	double r = 1.01 * r_g;
	ChVector<> hdims(hdimX - r, hdimY - r, 0);
	ChVector<> center(0, 0, 2 * r);

	while (gen.getTotalNumBodies() < num_particles) {
		gen.createObjectsBox(utils::POISSON_DISK, 2 * r, center, hdims);
		center.z += 2 * r;
	}

	std::cout << "Created " << gen.getTotalNumBodies() << " particles." << std::endl;

	return center.z;
}

// =============================================================================
// Utility function for displaying an ASCII progress bar for the quantity x
// which must be a value between 0 and n. The width 'w' represents the number
// of '=' characters corresponding to 100%.

void progressbar(unsigned int x, unsigned int n, unsigned int w = 50) {
	if ((x != n) && (x % (n / 100 + 1) != 0))
		return;

	float ratio = x / (float)n;
	unsigned int c = (unsigned int)(ratio * w);

	std::cout << std::setw(3) << (int)(ratio * 100) << "% [";
	for (unsigned int x = 0; x < c; x++)
		std::cout << "=";
	for (unsigned int x = c; x < w; x++)
		std::cout << " ";
	std::cout << "]\r" << std::flush;
}

// =============================================================================
int main(int argc, char* argv[]) {

	// --------------
	// Create system.
	// --------------
	// ----  Parallel
	std::cout << "Create Parallel DVI system" << std::endl;
	ChSystemParallelDVI* system = new ChSystemParallelDVI();

	system->Set_G_acc(ChVector<>(0, 0, -9.81));


	// ---------------------
	// Edit system settings.
	// ---------------------
	// Set solver parameters
	system->GetSettings()->solver.tolerance = tolerance;
	system->GetSettings()->solver.solver_mode = SLIDING;
	system->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
	system->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
	system->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
	system->GetSettings()->solver.max_iteration_bilateral = 1000;  // make 1000, should be about 220
	system->GetSettings()->solver.compute_N = false;
	system->GetSettings()->solver.alpha = 0;
	system->GetSettings()->solver.cache_step_length = true;
	system->GetSettings()->solver.use_full_inertia_tensor = false;
	system->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
	system->GetSettings()->solver.bilateral_clamp_speed = 1e8;
	system->ChangeSolverType(BB);
	system->SetLoggingLevel(LOG_INFO);
	system->SetLoggingLevel(LOG_TRACE);

	system->GetSettings()->collision.collision_envelope = 0.1 * r_g;


	system->GetSettings()->collision.bins_per_axis = vec3(100, 20, 25);
	system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;
	system->GetSettings()->collision.fixed_bins = true;

	// -------------------
	// Create the terrain.
	// -------------------

	CreateContainer(system);

	dof_container = new ChFluidContainer(system);
	CreateFluid(system);


	// --------------------------
	// Construct the M113 vehicle
	// --------------------------

	m_speedPIDVehicle.SetGains(4.0, 1.0, 0.00);

	// Create and initialize vehicle system
	M113_Vehicle vehicle(true, TrackShoeType::SINGLE_PIN, system);
	////vehicle.SetStepsize(0.0001);

	vehicle.Initialize(ChCoordsys<>(initLoc, initRot));

	vehicle.SetChassisVisualizationType(VisualizationType::MESH);
	vehicle.SetSprocketVisualizationType(VisualizationType::MESH);
	vehicle.SetIdlerVisualizationType(VisualizationType::MESH);
	vehicle.SetRoadWheelAssemblyVisualizationType(VisualizationType::MESH);
	vehicle.SetTrackShoeVisualizationType(VisualizationType::MESH);

	////vehicle.SetCollide(TrackCollide::NONE);
	////vehicle.SetCollide(TrackCollide::WHEELS_LEFT | TrackCollide::WHEELS_RIGHT);
	////vehicle.SetCollide(TrackCollide::ALL & (~TrackCollide::SPROCKET_LEFT) & (~TrackCollide::SPROCKET_RIGHT));

	// Create the powertrain system
	M113_SimplePowertrain powertrain;
	powertrain.Initialize(vehicle.GetChassisBody(), vehicle.GetDriveshaft());

	// ---------------
	// Simulation loop
	// ---------------

#ifdef CHRONO_OPENGL
	// Initialize OpenGL
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	gl_window.Initialize(1280, 720, "M113", system);
	gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
	gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

	// Number of simulation steps between two 3D view render frames
	int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

	// Run simulation for specified time.
	double time = 0;
	int sim_frame = 0;
	int out_frame = 0;
	int next_out_frame = 0;
	double exec_time = 0;
	int num_contacts = 0;

	// Inter-module communication data
	BodyStates shoe_states_left(vehicle.GetNumTrackShoes(LEFT));
	BodyStates shoe_states_right(vehicle.GetNumTrackShoes(RIGHT));
	TrackShoeForces shoe_forces_left(vehicle.GetNumTrackShoes(LEFT));
	TrackShoeForces shoe_forces_right(vehicle.GetNumTrackShoes(RIGHT));

	bool set_time = true;


	while (time < time_end) {
		// Collect output data from modules
		double throttle_input = 0;
		double steering_input = 0;
		double braking_input = 0;

		{

			double out_speed = m_speedPIDVehicle.Advance(vehicle, target_speed, time_step);
			if (time > .5) {
				ChClampValue(out_speed, -2.0, 2.0);
			}
			else {
				out_speed = 0;
			}

			if (out_speed > 0) {
				braking_input = 0;
				throttle_input = out_speed;
			}
			else {
				// Vehicle moving too fast: apply brakes
				braking_input = -out_speed;
				throttle_input = 0;
			}

			// Stop vehicle at time_brakes [s]
			if (vehicle.GetChassis()->GetPos().x > dist_end) {
				braking_input = 1;
				throttle_input = 0;
				if (set_time) {
					real time_end_temp = time + 1;
					time_end = Max(time_end_temp, time_end);
					set_time = false;
				}
			}

		}




		double powertrain_torque = powertrain.GetOutputTorque();
		double driveshaft_speed = vehicle.GetDriveshaftSpeed();
		vehicle.GetTrackShoeStates(LEFT, shoe_states_left);
		vehicle.GetTrackShoeStates(RIGHT, shoe_states_right);

		// Output
		if (sim_frame == next_out_frame) {
			cout << endl;
			cout << "---- Frame:          " << out_frame + 1 << endl;
			cout << "     Sim frame:      " << sim_frame << endl;
			cout << "     Time:           " << time << endl;
			cout << "     Avg. contacts:  " << num_contacts / out_steps << endl;
			cout << "     Throttle input: " << throttle_input << endl;
			cout << "     Braking input:  " << braking_input << endl;
			cout << "     Steering input: " << steering_input << endl;
			cout << "     Execution time: " << exec_time << endl;


			out_frame++;
			next_out_frame += out_steps;
			num_contacts = 0;
		}

		// Release the vehicle chassis at the end of the hold time.
		if (vehicle.GetChassisBody()->GetBodyFixed() && time > time_hold) {
			std::cout << std::endl << "Release vehicle t = " << time << std::endl;
			vehicle.GetChassisBody()->SetBodyFixed(false);
		}

		// Update modules (process inputs from other modules)
		powertrain.Synchronize(time, throttle_input, driveshaft_speed);
		vehicle.Synchronize(time, steering_input, braking_input, powertrain_torque, shoe_forces_left, shoe_forces_right);

		// Advance simulation for one timestep for all modules
		powertrain.Advance(time_step);
		vehicle.Advance(time_step);

#ifdef CHRONO_OPENGL
		if (gl_window.Active())
			gl_window.Render();
		else
			break;
#endif

		progressbar(out_steps + sim_frame - next_out_frame + 1, out_steps);

		// Periodically display maximum constraint violation
		if (monitor_bilaterals && sim_frame % bilateral_frame_interval == 0) {
			vehicle.LogConstraintViolations();
		}

		// Update counters.
		time += time_step;
		sim_frame++;
		exec_time += system->GetTimerStep();
		num_contacts += system->GetNcontacts();
	}

	// Final stats
	std::cout << "==================================" << std::endl;
	std::cout << "Simulation time:   " << exec_time << std::endl;
	std::cout << "Number of threads: " << threads << std::endl;

	return 0;
}
