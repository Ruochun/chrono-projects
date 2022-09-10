// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: RZ
// =============================================================================
//
// Demo to show Viper Rover operated on DEM Terrain
//
// =============================================================================

#include "chrono_models/robot/viper/Viper.h"

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_vehicle/ChVehicleModelData.h"

#include "chrono_thirdparty/filesystem/path.h"

#include <DEM/API.h>
//#include <core/ApiVersion.h>
//#include <core/utils/ThreadManager.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <map>
#include <random>
#include <cmath>

using namespace smug;
using namespace std::filesystem;

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::viper;

bool output = true;
const std::string out_dir = GetChronoOutputPath() + "SCM_DEF_SOIL";

// SCM grid spacing
double mesh_resolution = 0.02;

// Enable/disable bulldozing effects
bool enable_bulldozing = true;

// Enable/disable moving patch feature
bool enable_moving_patch = true;

// If true, use provided callback to change soil properties based on location
bool var_params = true;

// Define Viper rover wheel type
ViperWheelType wheel_type = ViperWheelType::RealWheel;

// Use custom material for the Viper Wheel
bool use_custom_mat = false;

// ChVector to/from float3
inline float3 ChVec2Float(const ChVector<>& vec) {
    return make_float3(vec.x(), vec.y(), vec.z());
}
inline ChVector<> Float2ChVec(float3 f3) {
    return ChVector<>(f3.x, f3.y, f3.z);
}

inline float4 ChQ2Float(const ChQuaternion<>& Q) {
    return make_float4(Q.e0(), Q.e1(), Q.e2(), Q.e3());
}

// Return customized wheel material parameters
std::shared_ptr<ChMaterialSurface> CustomWheelMaterial(ChContactMethod contact_method) {
    float mu = 0.4f;   // coefficient of friction
    float cr = 0.1f;   // coefficient of restitution
    float Y = 2e7f;    // Young's modulus
    float nu = 0.3f;   // Poisson ratio
    float kn = 2e5f;   // normal stiffness
    float gn = 40.0f;  // normal viscous damping
    float kt = 2e5f;   // tangential stiffness
    float gt = 20.0f;  // tangential viscous damping

    switch (contact_method) {
        case ChContactMethod::NSC: {
            auto matNSC = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            matNSC->SetFriction(mu);
            matNSC->SetRestitution(cr);
            return matNSC;
        }
        case ChContactMethod::SMC: {
            auto matSMC = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            matSMC->SetFriction(mu);
            matSMC->SetRestitution(cr);
            matSMC->SetYoungModulus(Y);
            matSMC->SetPoissonRatio(nu);
            matSMC->SetKn(kn);
            matSMC->SetGn(gn);
            matSMC->SetKt(kt);
            matSMC->SetGt(gt);
            return matSMC;
        }
        default:
            return std::shared_ptr<ChMaterialSurface>();
    }
}

int main(int argc, char* argv[]) {

    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    auto vis_mesh_file = GetChronoDataFile(".");

    std::cout << "Data dir is " << vis_mesh_file << std::endl;

    // Global parameter for moving patch size:
    double wheel_range = 0.5;
    ////double body_range = 1.2;

    // Create a Chrono::Engine physical system
    ChSystemSMC sys;
    ChVector<double> G = ChVector<double>(0, 0, -9.81);
    sys.Set_G_acc(G);

    const int nW = 4; // 4 wheels

    // Create the rover
    auto driver = chrono_types::make_shared<ViperDCMotorControl>();

    Viper viper(&sys, wheel_type);

    viper.SetDriver(driver);
    if (use_custom_mat)
        viper.SetWheelContactMaterial(CustomWheelMaterial(ChContactMethod::NSC));

    viper.Initialize(ChFrame<>(ChVector<>(-0.5, -0.0, -0.16), QUNIT));

    // Get wheels and bodies to set up SCM patches
    std::vector<std::shared_ptr<ChBodyAuxRef>> Wheels;
    std::vector<ChVector<>> wheel_pos;
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_LF)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_RF)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_LB)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_RB)->GetBody());
    
    auto Body_1 = viper.GetChassis()->GetBody();
    double total_mass = viper.GetRoverMass();
    double wheel_mass = viper.GetWheelMass();

    for (int i = 0; i < nW; i++) {
        wheel_pos.push_back(Wheels[i]->GetFrame_REF_to_abs().GetPos());
    }

    //////////////////////////////////////////////
    // Now step up SMUG
    //////////////////////////////////////////////

    DEMSolver DEM_sim;
    DEM_sim.SetVerbosity(INFO);
    DEM_sim.SetOutputFormat(DEM_OUTPUT_FORMAT::CSV);
    // DEM_sim.SetOutputContent(DEM_OUTPUT_CONTENT::FAMILY);
    DEM_sim.SetOutputContent(DEM_OUTPUT_CONTENT::XYZ);

    srand(759);

    float kg_g_conv = 1;
    // Define materials
    auto mat_type_terrain = DEM_sim.LoadMaterial({{"E", 2e9 * kg_g_conv}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.0}});
    auto mat_type_wheel = DEM_sim.LoadMaterial({{"E", 1e9 * kg_g_conv}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.0}});

    // Define the simulation world
    double world_y_size = 2.0;
    DEM_sim.InstructBoxDomainNumVoxel(22, 21, 21, (world_y_size) / std::pow(2, 16) / std::pow(2, 21));
    float bottom = -0.5;
    DEM_sim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);
    DEM_sim.AddBCPlane(make_float3(0, world_y_size / 2, 0), make_float3(0, -1, 0), mat_type_terrain);
    DEM_sim.AddBCPlane(make_float3(0, -world_y_size / 2, 0), make_float3(0, 1, 0), mat_type_terrain);
    // X-dir bounding planes
    DEM_sim.AddBCPlane(make_float3(-world_y_size * 2 / 2, 0, 0), make_float3(1, 0, 0), mat_type_terrain);
    DEM_sim.AddBCPlane(make_float3(world_y_size * 2 / 2, 0, 0), make_float3(-1, 0, 0), mat_type_terrain);

    // Define the wheel geometry
    float wheel_rad = 0.25;
    float wheel_width = 0.25;
    wheel_mass *= kg_g_conv;  // in kg or g
    // Our shelf wheel geometry is lying flat on ground with z being the axial direction
    float wheel_IYY = wheel_mass * wheel_rad * wheel_rad / 2;
    float wheel_IXX = (wheel_mass / 12) * (3 * wheel_rad * wheel_rad + wheel_width * wheel_width);
    float3 wheel_MOI = make_float3(wheel_IXX, wheel_IYY, wheel_IXX);
    auto wheel_template = DEM_sim.LoadClumpType(wheel_mass, wheel_MOI, "../data/clumps/ViperWheelSimple.csv", mat_type_wheel);
    // The file contains no wheel particles size info, so let's manually set them
    wheel_template->radii = std::vector<float>(wheel_template->nComp, 0.01);
    // This wheel template is `lying down', but our reported MOI info is assuming it's in a position to roll 
    // along X direction. Let's make it clear its principal axes is not what we used to report its component 
    // sphere relative positions.
    wheel_template->InformCentroidPrincipal(make_float3(0), make_float4(0.7071, 0, 0, 0.7071));

    // Then the ground particle template
    DEMClumpTemplate shape_template;
    shape_template.ReadComponentFromFile("../data/clumps/triangular_flat.csv");
    // Calculate its mass and MOI
    float mass = 2.6e3 * 5.5886717 * kg_g_conv;  // in kg or g
    float3 MOI = make_float3(1.8327927, 2.1580013, 0.77010059) * 2.6e3 * kg_g_conv;
    // Scale the template we just created
    std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates;
    std::vector<double> scales = {0.0014, 0.00063, 0.00033, 0.00022, 0.00015, 0.00009};
    std::for_each(scales.begin(), scales.end(), [](double& r) { r *= 10.; });
    for (double scaling : scales) {
        auto this_template = shape_template;
        this_template.mass = (double)mass * scaling * scaling * scaling;
        this_template.MOI.x = (double)MOI.x * (double)(scaling * scaling * scaling * scaling * scaling);
        this_template.MOI.y = (double)MOI.y * (double)(scaling * scaling * scaling * scaling * scaling);
        this_template.MOI.z = (double)MOI.z * (double)(scaling * scaling * scaling * scaling * scaling);
        std::cout << "Mass: " << this_template.mass << std::endl;
        std::cout << "MOIX: " << this_template.MOI.x << std::endl;
        std::cout << "MOIY: " << this_template.MOI.y << std::endl;
        std::cout << "MOIZ: " << this_template.MOI.z << std::endl;
        std::cout << "=====================" << std::endl;
        std::for_each(this_template.radii.begin(), this_template.radii.end(), [scaling](float& r) { r *= scaling; });
        std::for_each(this_template.relPos.begin(), this_template.relPos.end(), [scaling](float3& r) { r *= scaling; });
        this_template.materials = std::vector<std::shared_ptr<DEMMaterial>>(this_template.nComp, mat_type_terrain);
        ground_particle_templates.push_back(DEM_sim.LoadClumpType(this_template));
    }

    // Now we load part1 clump locations from a output file
    std::cout << "Making terrain..." << std::endl;
    auto part1_clump_xyz = DEM_sim.ReadClumpXyzFromCsv("GRC_10e6.csv");
    auto part1_clump_quaternion = DEM_sim.ReadClumpQuatFromCsv("GRC_10e6.csv");
    std::vector<float3> in_xyz;
    std::vector<float4> in_quat;
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
    unsigned int t_num = 0;
    for (int i = 0; i < scales.size(); i++) {
        // Our template names are 0001, 0002 etc.
        t_num++;
        char t_name[20];
        sprintf(t_name, "%04d", t_num);

        auto this_type_xyz = part1_clump_xyz[std::string(t_name)];
        auto this_type_quat = part1_clump_quaternion[std::string(t_name)];

        size_t n_clump_this_type = this_type_xyz.size();
        std::cout << "Loading clump " << std::string(t_name) << " which has particle num: " << n_clump_this_type << std::endl;
        // Prepare clump type identification vector for loading into the system (don't forget type 0 in
        // ground_particle_templates is the template for rover wheel)
        std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                 ground_particle_templates.at(t_num - 1));

        // Add them to the big long vector
        in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
        in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
        in_types.insert(in_types.end(), this_type.begin(), this_type.end());
        std::cout << "Added clump type " << t_num << std::endl;
    }
    // Finally, load the info into this batch
    DEMClumpBatch base_batch(in_xyz.size());
    base_batch.SetTypes(in_types);
    base_batch.SetPos(in_xyz);
    base_batch.SetOriQ(in_quat);

    DEM_sim.AddClumps(base_batch);

    /////////
    // Add wheel in DEM
    //////////
    
    // Instantiate this wheel
    std::cout << "Making wheels..." << std::endl;
    DEM_sim.SetFamilyFixed(100);
    std::vector<std::shared_ptr<DEMTracker>> trackers;
    std::vector<std::shared_ptr<DEMClumpBatch>> DEM_Wheels;
    for (int i = 0; i < nW; i++) {
        DEM_Wheels.push_back(DEM_sim.AddClumps(wheel_template, make_float3(wheel_pos[i].x(), wheel_pos[i].y(), wheel_pos[i].z())));
        // DEM_Wheel->SetOriQ(make_float4(0.7071, 0.7071, 0, 0));
        DEM_Wheels[i]->SetFamily(100);
        trackers.push_back(DEM_sim.Track(DEM_Wheels[i]));
    }

    //////
    // Make ready for DEM simulation
    ///////
    std::cout << "Begin initialization" << std::endl;
    auto max_v_finder = DEM_sim.CreateInspector("clump_max_absv");

    float base_step_size = 5e-7;
    float step_size = base_step_size;
    float base_vel = 0.4;
    DEM_sim.SetCoordSysOrigin("center");
    DEM_sim.SetInitTimeStep(step_size);
    DEM_sim.SetGravitationalAcceleration(ChVec2Float(G));
    // If you want to use a large UpdateFreq then you have to expand spheres to ensure safety
    DEM_sim.SetCDUpdateFreq(10);
    // DEM_sim.SetExpandFactor(1e-3);
    DEM_sim.SetMaxVelocity(15.0);
    DEM_sim.SetExpandSafetyParam(1.1);
    DEM_sim.SetInitBinSize(scales.at(2));
    
    DEM_sim.Initialize();
    for (const auto& tracker : trackers) {
        std::cout << "A tracker is tracking owner " << tracker->obj->ownerID << std::endl;
    }
    std::cout << "End initialization" << std::endl;

    ///////////////////////////////////////////
    // Real simulation
    ///////////////////////////////////////////

    float time_end = 2.0;
    unsigned int fps = 20;
    unsigned int report_freq = 50000;
    unsigned int param_update_freq = 50000;
    // unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    float frame_accu_thres = 1.0 / fps;
    unsigned int report_steps = (unsigned int)(1.0 / (report_freq * step_size));
    unsigned int param_update_steps = (unsigned int)(1.0 / (param_update_freq * step_size));

    path out_dir = current_path();
    out_dir += "/Viper_on_GRC";
    create_directory(out_dir);
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

    std::vector<ChQuaternion<>> wheel_rot(4);
    float max_v;
    int change_step = 0;
    float frame_accu = frame_accu_thres;
    for (float t = 0; t < time_end; t += step_size, curr_step++) {
        for (int i = 0; i < nW; i++) {
            wheel_pos[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos();
            trackers[i]->SetPos(ChVec2Float(wheel_pos[i]));
            wheel_rot[i] = Wheels[i]->GetFrame_REF_to_abs().GetRot();
            trackers[i]->SetOriQ(ChQ2Float(wheel_rot[i]));
            
            if (curr_step % report_steps == 0) {
                
                std::cout << "Wheel " << i << " position: " << wheel_pos[i].x() << ", " 
                      << wheel_pos[i].y() << ", " << wheel_pos[i].z() << std::endl;
                
                std::cout << "Wheel " << i << " rotation: " << wheel_rot[i].e0() << ", " 
                      << wheel_rot[i].e1() << ", " << wheel_rot[i].e2() << ", " << wheel_rot[i].e3() << std::endl;
            }
        }

        // if (curr_step % out_steps == 0) {
        if (frame_accu >= frame_accu_thres) {
            frame_accu = 0.;
            std::cout << "Frame: " << currframe << std::endl;
            DEM_sim.ShowThreadCollaborationStats();
            char filename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe++);
            DEM_sim.WriteSphereFile(std::string(filename));
            DEM_sim.ShowThreadCollaborationStats();
        }
        // Run DEM first
        DEM_sim.DoDynamics(step_size);

        // Then feed force
        for (int i = 0; i < nW; i++) {
            float3 F = trackers[i]->ContactAcc();
            F *= wheel_mass;
            float3 tor = trackers[i]->ContactAngAccLocal();
            tor = wheel_MOI * tor;
            Wheels[i]->Empty_forces_accumulators();
            Wheels[i]->Accumulate_force(Float2ChVec(F), wheel_pos[i], false);
            Wheels[i]->Accumulate_torque(Float2ChVec(tor), true); // torque in SMUG is local
        }
        sys.DoStepDynamics(step_size);
        viper.Update();
        t += step_size;
        frame_accu += step_size;

        // if (curr_step % param_update_steps == 0 && t < 0.2) {
        // // if (t < 0.25) {
        //     DEM_sim.DoDynamicsThenSync(0);
        //     max_v = max_v_finder->GetValue();
        //     float multiplier = max_v / base_vel;
        //     step_size = base_step_size / multiplier;
        //     DEM_sim.SetInitTimeStep(step_size);
        //     // DEM_sim.SetMaxVelocity(max_v * 1.2);
        //     DEM_sim.UpdateSimParams();
        //     std::cout << "Max vel in simulation is " << max_v << std::endl;
        //     std::cout << "Step size in simulation is " << step_size << std::endl;
        // }
        if (t > 0.2 && change_step == 0) {
            DEM_sim.DoDynamicsThenSync(0);
            step_size = 1e-6;
            DEM_sim.SetInitTimeStep(step_size);
            DEM_sim.UpdateSimParams();
            change_step = 1;
        } else if (t > 0.3 && change_step == 1) {
            DEM_sim.DoDynamicsThenSync(0);
            step_size = 2e-6;
            DEM_sim.SetInitTimeStep(step_size);
            DEM_sim.UpdateSimParams();
            change_step = 2;
        } else if (t > 0.4 && change_step == 2) {
            DEM_sim.DoDynamicsThenSync(0);
            step_size = 5e-6;
            DEM_sim.SetInitTimeStep(step_size);
            DEM_sim.UpdateSimParams();
            change_step = 3;
        }


        if (curr_step % report_steps == 0) {
            float3 body_pos = ChVec2Float(Body_1->GetFrame_REF_to_abs().GetPos());
            std::cout << "Rover body is at " << body_pos.x << ", " << body_pos.y << ", " << body_pos.z << std::endl;
            std::cout << "Time is " << t << std::endl;
            max_v = max_v_finder->GetValue();
            std::cout << "Max vel in simulation is " << max_v << std::endl;
            std::cout << "========================" << std::endl;
        }
    }

    DEM_sim.ShowThreadCollaborationStats();
    DEM_sim.ShowTimingStats();

    return 0;
}
