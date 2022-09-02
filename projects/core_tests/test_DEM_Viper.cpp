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
// Demo to show Viper Rover operated on SCM Terrain
//
// =============================================================================

#include "chrono_models/robot/viper/Viper.h"

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_vehicle/ChVehicleModelData.h"

#include "chrono_thirdparty/filesystem/path.h"

#include <DEM/API.h>
#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <map>
#include <random>
#include <cmath>

using namespace sgps;
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

    // Create the rover
    auto driver = chrono_types::make_shared<ViperDCMotorControl>();

    Viper viper(&sys, wheel_type);

    viper.SetDriver(driver);
    if (use_custom_mat)
        viper.SetWheelContactMaterial(CustomWheelMaterial(ChContactMethod::NSC));

    viper.Initialize(ChFrame<>(ChVector<>(-5, 0, -0.2), QUNIT));

    // Get wheels and bodies to set up SCM patches
    auto Wheel_1 = viper.GetWheel(ViperWheelID::V_LF)->GetBody();
    auto Wheel_2 = viper.GetWheel(ViperWheelID::V_RF)->GetBody();
    auto Wheel_3 = viper.GetWheel(ViperWheelID::V_LB)->GetBody();
    auto Wheel_4 = viper.GetWheel(ViperWheelID::V_RB)->GetBody();
    auto Body_1 = viper.GetChassis()->GetBody();
    double total_mass = viper.GetRoverMass();

    //////////////////////////////////////////////
    // Now step up SGPS
    //////////////////////////////////////////////

    DEMSolver DEM_sim;
    DEM_sim.SetVerbosity(INFO);
    DEM_sim.SetOutputFormat(DEM_OUTPUT_FORMAT::CSV);
    // DEM_sim.SetOutputContent(DEM_OUTPUT_CONTENT::FAMILY);
    DEM_sim.SetOutputContent(DEM_OUTPUT_CONTENT::XYZ);

    srand(759);

    float kg_g_conv = 1;
    // Define materials
    auto mat_type_terrain = DEM_sim.LoadMaterialType(2e9 * kg_g_conv, 0.3, 0.3, 0.5, 0.0);
    auto mat_type_wheel = DEM_sim.LoadMaterialType(1e9 * kg_g_conv, 0.3, 0.3, 0.5, 0.0);

    // Define the simulation world
    double world_y_size = 0.99;
    DEM_sim.InstructBoxDomainNumVoxel(21, 21, 22, world_y_size / std::pow(2, 16) / std::pow(2, 21));
    // Add 5 bounding planes around the simulation world, and leave the top open
    DEM_sim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);
    float bottom = -0.5;
    DEM_sim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);

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
    auto part1_clump_xyz = DEM_sim.ReadClumpXyzFromCsv("GRC_1e6.csv");
    auto part1_clump_quaternion = DEM_sim.ReadClumpQuatFromCsv("GRC_1e6.csv");
    std::vector<float3> in_xyz;
    std::vector<float4> in_quat;
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
    unsigned int t_num;
    for (int i = 0; i < scales.size(); i++) {
        // Our template names are 0001, 0002 etc.
        t_num++;
        char t_name[20];
        sprintf(t_name, "%04d", t_num);

        auto this_type_xyz = part1_clump_xyz[std::string(t_name)];
        auto this_type_quat = part1_clump_quaternion[std::string(t_name)];

        size_t n_clump_this_type = this_type_xyz.size();
        // Prepare clump type identification vector for loading into the system (don't forget type 0 in
        // ground_particle_templates is the template for rover wheel)
        std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                 ground_particle_templates.at(t_num - 1));

        // Add them to the big long vector
        in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
        in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
        in_types.insert(in_types.end(), this_type.begin(), this_type.end());
    }
    // Finally, load the info into this batch
    DEMClumpBatch base_batch(in_xyz.size());
    base_batch.SetTypes(in_types);
    base_batch.SetPos(in_xyz);
    base_batch.SetOriQ(in_quat);

    DEM_sim.AddClumps(base_batch);
    

    // Make ready for simulation
    float step_size = 5e-7;
    DEM_sim.SetCoordSysOrigin("center");
    DEM_sim.SetInitTimeStep(step_size);
    DEM_sim.SetGravitationalAcceleration(make_float3(0, 0, -9.8));
    // If you want to use a large UpdateFreq then you have to expand spheres to ensure safety
    DEM_sim.SetCDUpdateFreq(10);
    // DEM_sim.SetExpandFactor(1e-3);
    DEM_sim.SetMaxVelocity(20.);
    DEM_sim.SetExpandSafetyParam(1.2);
    DEM_sim.SetInitBinSize(scales.at(2));
    DEM_sim.Initialize();

    ///////////////////////////////////////////
    // Real simulation
    ///////////////////////////////////////////

    float t = 0.0;
    while (t < 1.0) {

        if (output) {
            // std::cout << sys.GetChTime() << viper.GetWheelTracTorque(ViperWheelID::V_LF)
            //     << viper.GetWheelTracTorque(ViperWheelID::V_RF) << viper.GetWheelTracTorque(ViperWheelID::V_LB)
            //     << viper.GetWheelTracTorque(ViperWheelID::V_RB) << std::endl;
            ChFrame<> body_ref_frame = Wheel_1->GetFrame_REF_to_abs();
            ChVector<> body_pos = body_ref_frame.GetPos();  
            std::cout << "Wheel pos: " << body_pos.x() << " " << body_pos.y() << " " << body_pos.z() << std::endl;
        }

        ChVector<double> quarter_force = G * (-total_mass / 4);
        Wheel_1->Empty_forces_accumulators();
        Wheel_1->Accumulate_force(quarter_force, Wheel_1->GetFrame_REF_to_abs().GetPos(), false);
        // Wheel_1->Accumulate_torque(Body_torque, false);
        Wheel_2->Empty_forces_accumulators();
        Wheel_2->Accumulate_force(quarter_force, Wheel_2->GetFrame_REF_to_abs().GetPos(), false);
        // Wheel_2->Accumulate_torque(Body_torque, false);
        Wheel_3->Empty_forces_accumulators();
        Wheel_3->Accumulate_force(quarter_force, Wheel_3->GetFrame_REF_to_abs().GetPos(), false);
        // Wheel_3->Accumulate_torque(Body_torque, false);
        Wheel_4->Empty_forces_accumulators();
        Wheel_4->Accumulate_force(quarter_force, Wheel_4->GetFrame_REF_to_abs().GetPos(), false);
        // Wheel_4->Accumulate_torque(Body_torque, false);

        sys.DoStepDynamics(5e-4);
        viper.Update();
        t += 5e-4;
        ////terrain.PrintStepStatistics(std::cout);
    }

    return 0;
}
