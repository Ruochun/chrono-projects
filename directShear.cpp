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
// Author: Daniel Melanz, Radu Serban
// =============================================================================
//
// ChronoParallel demo program for shearing studies.
//
// The system contains a shearing box composed of three bodies: (1) a containing
// bin, (2) a shearing plate, and (3) a load body. Granular material sits inside
// of the containing bin and shearing plate and is compressed from the top by
// the load body. During the shearing mode, the shear plate is translated in the
// x-direction at a specified velocity.
//
// The global reference frame has Z up.
// All units SI (CGS, i.e., centimeter - gram - second)
//
// =============================================================================

#include <stdio.h>
#include <vector>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "chrono_parallel/ChSystemParallel.h"
#include "chrono_parallel/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

//#undef CHRONO_PARALLEL_HAS_OPENGL
#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Problem definitions
// -----------------------------------------------------------------------------

// Comment the following line to use DVI contact
#define DEM

enum ProblemType {
  SETTLING,
  PRESSING,
  SHEARING,
  TESTING
};

ProblemType problem = TESTING;

// -----------------------------------------------------------------------------
// Conversion factors
// -----------------------------------------------------------------------------

// Conversion for Young's modulus and pressure
// [Y] = Pa = N / m^2 = kg / m / s^2
double Pa2cgs = 10;

// Conversion for dissipation factor alpha
// [alpha] = s / m
double alpha2cgs = 0.01;
 
// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Perform shearing action via a linear actuator or kinematically?
bool use_actuator = true;

// Simulation parameters
double gravity = 981;  // [cm/s^2] 

double time_settling_min = 0.1;
double time_settling_max = 0.2;

double time_pressing_min = 0.1;
double time_pressing_max = 0.2;

double time_shearing = 10;  //0.06;

#ifdef DEM
double time_step = 1e-5;
#else
double time_step = 1e-4;
int max_iteration_normal = 50;
int max_iteration_sliding = 100;
int max_iteration_spinning = 0;
float contact_recovery_speed = 1;
#endif

int max_iteration_bilateral = 100;
double tolerance = 1e-4;

// Output
#ifdef DEM
const std::string out_dir = "../DIRECTSHEAR_DEM/";
#else
const std::string out_dir = "../DIRECTSHEAR_DVI/";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string shear_file = out_dir + "/shear.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string settled_ckpnt_file = out_dir + "/settled.dat";
const std::string pressed_ckpnt_file = out_dir + "/pressed.dat";

int out_fps_settling = 120;
int out_fps_pressing = 120;
int out_fps_shearing = 120;

int timing_frame = 10;   // output detailed step timing at this frame

// Parameters for the containing bin and shear box
int        Id_ground = -1;             // body ID for the ground (containing bin)
int        Id_box = -2;                // body ID for the shear box

double     hdimX = 6.0 / 2;            // [cm] half-length in x direction
double     hdimY = 6.0 / 2;            // [cm] half-depth in y direction
double     hdimZ = 3.0 / 2;            // [cm] half-height in z direction
double     hthick = 1.0 / 2;           // [cm] half-thickness of the walls

float      Y_walls = Pa2cgs * 2e6;
float      alpha_walls = alpha2cgs * 0.4;
float      mu_walls = 0.3f;

// Parameters for the load plate
int        Id_plate = -3;              // body ID for the load plate

float      Y_plate = Pa2cgs * 2e6;
float      alpha_plate = alpha2cgs * 0.4;
float      mu_plate = 0.3f;

// Applied normal pressure
double     normalPressure = Pa2cgs * 16888.1;  // 16,888.1 Pa // 44,127.0 Pa// 71,365.9 Pa

// Desired shearing velocity [cm/s]
double     desiredVelocity = 0.5;  // 0.066;

// Parameters for the granular material
int        Id_g = 1;
double     r_g = 0.3;                    // [cm] radius of granular sphers
double     rho_g = 2.500;                // [g/cm^3] density of granules

double     desiredBulkDensity = 1.3894;  // [g/cm^3] desired bulk density

float      Y_g = Pa2cgs * 2e6;
float      alpha_g = alpha2cgs * 0.4;
float      mu_g = 0.5f;

// Parameters of the testing ball
int        Id_ball = -4;
double     mass_ball = 200;               // [g] mass of testing ball
double     radius_ball = 0.9 * hdimX;     // [cm] radius of testing ball


// =============================================================================
// Create the containing bin (the ground) and shear box.
// =============================================================================

ChSharedPtr<ChBody> CreateContainerAndBox(ChSystemParallel* system)
{
  // -------------------------------
  // Create a material for the walls
  // -------------------------------

#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_walls;
  mat_walls = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_walls->SetYoungModulus(Y_walls);
  mat_walls->SetFriction(mu_walls);
  mat_walls->SetDissipationFactor(alpha_walls);
#else
  ChSharedPtr<ChMaterialSurface> mat_walls(new ChMaterialSurface);
  mat_walls->SetFriction(mu_walls);
#endif

  // ----------------------
  // Create the ground body -- always FIRST body in system
  // ----------------------

#ifdef DEM
  ChSharedPtr<ChBodyDEM> ground(new ChBodyDEM(new ChCollisionModelParallel));
  ground->SetMaterialSurfaceDEM(mat_walls);
#else
  ChSharedPtr<ChBody> ground(new ChBody(new ChCollisionModelParallel));
  ground->SetMaterialSurface(mat_walls);
#endif

  ground->SetIdentifier(Id_ground);
  ground->SetBodyFixed(true);
  ground->SetCollide(true);

  // Attach geometry of the containing bin
  ground->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(ground.get_ptr(), ChVector<>(hdimX, hdimY, hthick), ChVector<>(0, 0, -hthick));
  utils::AddBoxGeometry(ground.get_ptr(), ChVector<>(hthick, hdimY, hdimZ), ChVector<>(-hdimX - hthick, 0, hdimZ));
  utils::AddBoxGeometry(ground.get_ptr(), ChVector<>(hthick, hdimY, hdimZ), ChVector<>( hdimX + hthick, 0, hdimZ));
  utils::AddBoxGeometry(ground.get_ptr(), ChVector<>(hdimX, hthick, hdimZ), ChVector<>(0, -hdimY - hthick, hdimZ));
  utils::AddBoxGeometry(ground.get_ptr(), ChVector<>(hdimX, hthick, hdimZ), ChVector<>(0,  hdimY + hthick, hdimZ));
  ground->GetCollisionModel()->BuildModel();

  system->AddBody(ground);

  // --------------------
  // Create the shear box -- always SECOND body in system
  // --------------------

  // Initially, the shear box is fixed to ground.
  // During the shearing phase it may be released (if using an actuator)
#ifdef DEM
  ChSharedBodyDEMPtr box(new ChBodyDEM(new ChCollisionModelParallel));
  box->SetMaterialSurfaceDEM(mat_walls);
#else
  ChSharedBodyPtr box(new ChBody(new ChCollisionModelParallel));
  box->SetMaterialSurface(mat_walls);
#endif

  box->SetIdentifier(Id_box);
  box->SetPos(ChVector<>(0, 0, 2 * hdimZ));
  box->SetCollide(true);
  box->SetBodyFixed(true);

  box->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(box.get_ptr(), ChVector<>(hthick, hdimY, 5 * hdimZ), ChVector<>(-hdimX - hthick, 0, 5 * hdimZ));
  utils::AddBoxGeometry(box.get_ptr(), ChVector<>(hthick, hdimY, 5 * hdimZ), ChVector<>( hdimX + hthick, 0, 5 * hdimZ));
  utils::AddBoxGeometry(box.get_ptr(), ChVector<>(hdimX, hthick, 5 * hdimZ), ChVector<>(0, -hdimY - hthick, 5 * hdimZ));
  utils::AddBoxGeometry(box.get_ptr(), ChVector<>(hdimX, hthick, 5 * hdimZ), ChVector<>(0,  hdimY + hthick, 5 * hdimZ));
  box->GetCollisionModel()->BuildModel();

  system->AddBody(box);

  // -------------------------------------------------------------------------
  // Create a translational joint between ground and box and a linear actuator
  // -------------------------------------------------------------------------

  if (use_actuator) {
    ChSharedPtr<ChLinkLockPrismatic> prismatic(new ChLinkLockPrismatic);
    prismatic->Initialize(ground, box, ChCoordsys<>(ChVector<>(0, 0, 2 * hdimZ), Q_from_AngY(CH_C_PI_2)));
    system->AddLink(prismatic);

    ChSharedPtr<ChLinkLinActuator> actuator(new ChLinkLinActuator);
    ChVector<> pt1(0, 0, 2 * hdimZ);
    ChVector<> pt2(1, 0, 2 * hdimZ);
    actuator->SetName("actuator");
    actuator->Initialize(ground, box, false, ChCoordsys<>(pt1, QUNIT), ChCoordsys<>(pt2, QUNIT));
    actuator->Set_lin_offset(1);
    system->AddLink(actuator);
  }


  // Return a pointer to the shear box body
  return box;
}

// =============================================================================
// Create the load plate
// =============================================================================

ChSharedPtr<ChBody> CreateLoadPlate(ChSystemParallel* system, double z)
{
  // ------------------------------
  // Create a material for the body
  // ------------------------------

#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_plate;
  mat_plate = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_plate->SetYoungModulus(Y_plate);
  mat_plate->SetFriction(mu_plate);
  mat_plate->SetDissipationFactor(alpha_plate);
#else
  ChSharedPtr<ChMaterialSurface> mat_plate(new ChMaterialSurface);
  mat_plate->SetFriction(mu_plate);
#endif

  // ---------------------
  // Create the plate body
  // ---------------------

  // Set plate dimension, shrinking by a half-radius on each side
  double hdimX_p = hdimX - r_g / 2;
  double hdimY_p = hdimY - r_g / 2;
  double hdimZ_p = hdimZ;

  // Estimate plate mass from desired applied normal pressure
  double area = 4 * hdimX_p * hdimY_p;
  double mass = normalPressure * area / gravity;

#ifdef DEM
  ChSharedBodyDEMPtr plate(new ChBodyDEM(new ChCollisionModelParallel));
  plate->SetMaterialSurfaceDEM(mat_plate);
#else
  ChSharedBodyPtr plate(new ChBody(new ChCollisionModelParallel));
  plate->SetMaterialSurface(mat_plate);
#endif

  plate->SetIdentifier(Id_plate);
  plate->SetMass(mass);
  plate->SetPos(ChVector<>(0, 0, z + hdimZ_p));
  plate->SetCollide(true);
  plate->SetBodyFixed(false);

  plate->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(plate.get_ptr(), ChVector<>(hdimX_p, hdimY_p, hdimZ_p));
  plate->GetCollisionModel()->BuildModel();

  system->AddBody(plate);

  // --------------------------------------------------------
  // Create a translational joint between shear box and plate
  // --------------------------------------------------------

  ChSharedPtr<ChBody> box(system->Get_bodylist()->at(1));

  ChSharedPtr<ChLinkLockPrismatic> prismatic(new ChLinkLockPrismatic);
  prismatic->Initialize(box, plate, ChCoordsys<>(ChVector<>(0, 0, z + hdimZ_p), QUNIT));
  system->AddLink(prismatic);

  // Return a pointer to the plate body
  return plate;
}


// =============================================================================
// Create the granular material
//
// Granular material consisting of identical spheres with specified radius and
// material properties; the spheres are generated in a number of vertical
// layers with locations within each layer obtained using Poisson Disk sampling,
// thus ensuring that no two spheres are closer than twice the radius.
// =============================================================================

int CreateObjects(ChSystemParallel* system)
{
  // -------------------------------------------
  // Create a material for the granular material
  // -------------------------------------------

#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_g;
  mat_g = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_g->SetYoungModulus(Y_g);
  mat_g->SetFriction(mu_g);
  mat_g->SetDissipationFactor(alpha_g);
#else
  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
  mat_g->SetFriction(mu_g);
#endif

  // ---------------------------------------------
  // Create a mixture entirely made out of spheres
  // ---------------------------------------------

  // Create the particle generator with a mixture of 100% spheres
  utils::Generator gen(system);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
#ifdef DEM
  m1->setDefaultMaterialDEM(mat_g);
#else
  m1->setDefaultMaterialDVI(mat_g);
#endif
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(r_g);

  // Ensure that all generated particle bodies will have positive IDs.
  gen.setBodyIdentifier(Id_g);

  // ----------------------
  // Generate the particles
  // ----------------------

  double r = 1.01 * r_g;
  ChVector<> hdims(hdimX - r, hdimY - r, 0);
  ChVector<> center(0, 0, 2 * r);

  while (center.z < 6 * hdimZ)
  {
    gen.createObjectsBox(utils::POISSON_DISK, 2 * r, center, hdims);
    center.z += 2 * r;
  }

  // Return the number of generated particles.
  return gen.getTotalNumBodies();
}


// =============================================================================
// Create a single large sphere (for use in TESTING)
// =============================================================================
ChSharedPtr<ChBody> CreateBall(ChSystemParallel* system)
{
  // ------------------------------
  // Create a material for the ball
  // ------------------------------

#ifdef DEM
  ChSharedPtr<ChMaterialSurfaceDEM> mat_g;
  mat_g = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  mat_g->SetYoungModulus(Y_g);
  mat_g->SetFriction(mu_g);
  mat_g->SetDissipationFactor(alpha_g);
#else
  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
  mat_g->SetFriction(mu_g);
#endif

  // ---------------
  // Create the ball
  // ---------------

#ifdef DEM
  ChSharedBodyDEMPtr ball(new ChBodyDEM(new ChCollisionModelParallel));
  ball->SetMaterialSurfaceDEM(mat_g);
#else
  ChSharedBodyPtr ball(new ChBody(new ChCollisionModelParallel));
  ball->SetMaterialSurface(mat_g);
#endif

  ball->SetIdentifier(Id_ball);
  ball->SetMass(mass_ball);
  ball->SetPos(ChVector<>(0, 0, 1.01 * radius_ball));
  ball->SetCollide(true);
  ball->SetBodyFixed(false);

  ball->GetCollisionModel()->ClearModel();
  utils::AddSphereGeometry(ball.get_ptr(), radius_ball);
  ball->GetCollisionModel()->BuildModel();

  system->AddBody(ball);

  // Return a pointer to the ball body
  return ball;
}


// =============================================================================
// Find the height of the highest and lowest sphere in the granular mix. We only
// look at bodies whith positive identifiers (to exclude all other bodies).
// =============================================================================

void FindHeightRange(ChSystemParallel* sys, double& lowest, double& highest)
{
  highest = -1000;
  lowest = 1000;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() <= 0)
      continue;
    double h = body->GetPos().z;
    if (h < lowest)       lowest = h;
    else if (h > highest) highest = h;
  }
}


// =============================================================================
//
//// TODO:  cannot do this with DEM!!!!!
//
// =============================================================================

int setBulkDensity(ChSystem* sys, double bulkDensity)
{
  double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;

  double normalPlateHeight = sys->Get_bodylist()->at(1)->GetPos().z - hdimZ;
  double bottomHeight = 0;
  int numBodies = sys->Get_bodylist()->size();
  double boxVolume = hdimX * 2 * hdimX * 2 * (normalPlateHeight - bottomHeight);
  double granularVolume = (numBodies - 3)*vol_g;
  double reqDensity = bulkDensity*boxVolume / granularVolume;
  for (int i = 0; i < sys->Get_bodylist()->size(); ++i) {
    ChBody* body = (ChBody*)sys->Get_bodylist()->at(i);
    if (body->GetIdentifier() > 1) {
      body->SetMass(reqDensity*vol_g);
    }
  }

  cout << "N Bodies: " << numBodies << endl;
  cout << "Box Volume: " << boxVolume << endl;
  cout << "Granular Volume: " << granularVolume << endl;
  cout << "Desired bulk density = " << bulkDensity << ", Required Body Density = " << reqDensity << endl;

  return 0;
}


// =============================================================================

int main(int argc, char* argv[])
{
  // Create output directories.
  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if(ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

  // -------------
  // Create system
  // -------------

#ifdef DEM
  cout << "Create DEM system" << endl;
  ChSystemParallelDEM* msystem = new ChSystemParallelDEM();
#else
  cout << "Create DVI system" << endl;
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
#endif

  msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

  // Set number of threads.
  int max_threads = msystem->GetParallelThreadNumber();
  if (threads > max_threads) threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  cout << "Using " << threads << " threads" << endl;

  msystem->GetSettings()->max_threads = threads;
  msystem->GetSettings()->perform_thread_tuning = thread_tuning;

  // Edit system settings
  msystem->SetTol(tolerance);
  msystem->SetTolSpeeds(tolerance);
  msystem->SetStep(time_step);

  msystem->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
  msystem->GetSettings()->solver.tolerance = tolerance;

#ifdef DEM
  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_R;
#else
  msystem->GetSettings()->solver.solver_mode = SLIDING;
  msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
  msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
  msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
  msystem->GetSettings()->solver.alpha = 0;
  msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
  msystem->SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  msystem->ChangeSolverType(APGDBLAZE);

  msystem->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);
  msystem->GetSettings()->collision.min_body_per_bin = 50;
  msystem->GetSettings()->collision.max_body_per_bin = 100;

  // --------------
  // Problem set up
  // --------------

  // Depending on problem type:
  // - Select end simulation time
  // - Select output FPS
  // - Create / load objects

  double time_end;
  int out_fps;
  ChSharedPtr<ChBody> shearBox;
  ChSharedPtr<ChBody> loadPlate;
  ChSharedPtr<ChLinkLinActuator> actuator;

  switch (problem) {
  case SETTLING:
  {
    time_end = time_settling_max;
    out_fps = out_fps_settling;

    shearBox = CreateContainerAndBox(msystem);

    // Create granular material and containing bin
    cout << "Create granular material...  ";
    int num_particles = CreateObjects(msystem);
    cout << num_particles << " particles" << endl;

    break;
  }

  case PRESSING:
  {
    time_end = time_pressing_max;
    out_fps = out_fps_pressing;

    // Create bodies (container and shear box) from checkpoint file.
    cout << "Read checkpoint data from " << settled_ckpnt_file;
    utils::ReadCheckpoint(msystem, settled_ckpnt_file);
    cout << "  done.  Read " << msystem->Get_bodylist()->size() << " bodies." << endl;

    // Create the load plate just above the granular material.
    double highest, lowest;
    FindHeightRange(msystem, lowest, highest);
    loadPlate = CreateLoadPlate(msystem, highest + 2 * r_g);

    break;
  }

  case SHEARING:
  {
    time_end = time_shearing;
    out_fps = out_fps_shearing;

    // Create bodies (container, shear box, and plate) from checkpoint file.
    cout << "Read checkpoint data from " << pressed_ckpnt_file;
    utils::ReadCheckpoint(msystem, pressed_ckpnt_file);
    cout << "  done.  Read " << msystem->Get_bodylist()->size() << " bodies." << endl;

    // Fix the loadPlate and release the shearBox (if using an actuator)
    loadPlate = ChSharedPtr<ChBody>(msystem->Get_bodylist()->at(0));
    shearBox = ChSharedPtr<ChBody>(msystem->Get_bodylist()->at(1));
    loadPlate->SetBodyFixed(true);
    shearBox->SetBodyFixed(!use_actuator);

    // Grab the pointer to the actuator
    if (use_actuator)
      actuator = msystem->SearchLink("actuator").StaticCastTo<ChLinkLinActuator>();

    //setBulkDensity(msystem, desiredBulkDensity);

    break;
  }

  case TESTING:
  {
    time_end = 2;
    out_fps = 60;

    // Create container and shear box.
    shearBox = CreateContainerAndBox(msystem);

    // Create the test ball
    CreateBall(msystem);

    // Create the load plate just above the ball
    loadPlate = CreateLoadPlate(msystem, 2.1 * radius_ball);

    // Release the shear box (if using an actuator)
    shearBox->SetBodyFixed(!use_actuator);

    // Grab the pointer to the actuator
    if (use_actuator)
      actuator = msystem->SearchLink("actuator").StaticCastTo<ChLinkLinActuator>();

    break;
  }

  }

  // ----------------------
  // Perform the simulation
  // ----------------------

  // Set number of simulation steps and output frequency
  int num_steps = std::ceil(time_end / time_step);
  int out_steps = std::ceil((1.0 / time_step) / out_fps);

  // Initialize counters and output files.
  double time = 0;
  int sim_frame = 0;
  int out_frame = 0;
  int next_out_frame = 0;
  double exec_time = 0;
  int num_contacts = 0;
  ChStreamOutAsciiFile sfile(stats_file.c_str());
  ChStreamOutAsciiFile shearStream(shear_file.c_str());

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow &gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "Direct Shear Test", msystem);
  gl_window.SetCamera(ChVector<>(0,-10,0), ChVector<>(0,0,0),ChVector<>(0,0,1));
#endif

  // Loop until reaching the end time...
  while (time < time_end) {

    // Current position and velocity of the shear box
    const ChVector<>& pos_old = shearBox->GetPos();
    const ChVector<>& vel_old = shearBox->GetPos_dt();

    // Calculate running averages for particle min/max heights
    //// TODO
    double highest, lowest;
    FindHeightRange(msystem, lowest, highest);

    if (sim_frame == next_out_frame) {
      char filename[100];
      sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
      utils::WriteShapesPovray(msystem, filename, false);

      cout << "------------ Output frame:   " << out_frame << endl;
      cout << "             Sim frame:      " << sim_frame << endl;
      cout << "             Time:           " << time << endl;
      cout << "             Shear box pos:  " << pos_old.x << endl;
      cout << "             Lowest point:   " << lowest << endl;
      cout << "             Highest point:  " << highest << endl;
      cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
      cout << "             Execution time: " << exec_time << endl;

      sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";

      // Create a checkpoint from the current state.
      if (problem == SETTLING || problem == PRESSING) {
        cout << "             Write checkpoint data " << flush;
        if (problem == SETTLING) utils::WriteCheckpoint(msystem, settled_ckpnt_file);
        else                     utils::WriteCheckpoint(msystem, pressed_ckpnt_file);
        cout << msystem->Get_bodylist()->size() << " bodies" << endl;
      }

      // Increment counters
      out_frame++;
      next_out_frame += out_steps;
      num_contacts = 0;
    }

    ////if (problem == SETTLING && time > time_settling_min && CheckSettled(msystem, zero_v)) {
    ////  cout << "Granular material settled...  time = " << time << endl;
    ////  break;
    ////}

    // Advance simulation time
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(time_step);
      gl_window.Render();
    }
#else
    msystem->DoStepDynamics(time_step);
#endif

    if (problem == SHEARING || problem == TESTING) {
      // Calculate desired position of shear box at new time
      double xpos_new = pos_old.x + desiredVelocity * time_step;

      // Get/estimate the current reaction force and impose shear box position
      double cnstr_force = 0;
      if (use_actuator) {
        cnstr_force = actuator->Get_react_force().x;
        if (ChSharedPtr<ChFunction_Const> fun = actuator->Get_dist_funct().DynamicCastTo<ChFunction_Const>())
          fun->Set_yconst(xpos_new);
      }
      else {
        cnstr_force = (shearBox->GetPos_dt().x - vel_old.x) / time_step;
        shearBox->SetPos(ChVector<>(xpos_new, pos_old.y, pos_old.z));
        shearBox->SetPos_dt(ChVector<>(desiredVelocity, 0, 0));
      }

      ////shearStream << time << ", " << xpos_new << ", " << cnstr_force << ", \n";
      ////cout << "X pos: " << xpos_new << " X react: " << cnstr_force << endl;
    }

    // Increment counters
    time += time_step;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
    num_contacts += msystem->GetNcontacts();

    // If requested, output detailed timing information for this step
    if (sim_frame == timing_frame)
      msystem->PrintStepStats();
  }

  // ----------------
  // Final processing
  // ----------------

  // Create a checkpoint from the last state
  if (problem == SETTLING || problem == PRESSING) {
    cout << "             Write checkpoint data " << flush;
    if (problem == SETTLING) utils::WriteCheckpoint(msystem, settled_ckpnt_file);
    else                     utils::WriteCheckpoint(msystem, pressed_ckpnt_file);
    cout << msystem->Get_bodylist()->size() << " bodies" << endl;
  }

  // Final stats
  cout << "==================================" << endl;
  cout << "Number of bodies:  " << msystem->Get_bodylist()->size() << endl;
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}

