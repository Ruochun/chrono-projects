// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2019 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Nic Olsen
// =============================================================================
// Chrono::Granular simulation in which a cylinder is filled with granular
// material and then raised slightly, allowing material to flow out and around
// the cylinder for comparison with the analytical hydrostatic result.
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_granular/ChGranularData.h"
#include "chrono_granular/api/ChApiGranularChrono.h"
#include "chrono_granular/physics/ChGranular.h"
#include "chrono_granular/physics/ChGranularTriMesh.h"
#include "chrono_granular/utils/ChGranularJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::granular;

void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> <output_dir> <radius> <density>" << std::endl;
}

void writeMeshFrames(std::ostringstream& outstream,
                     const std::string obj_name,
                     const ChVector<>& pos,
                     const float3 mesh_scaling) {
    outstream << obj_name << ",";

    // Get basis vectors
    ChVector<> vx(1, 0, 0);
    ChVector<> vy(0, 1, 0);
    ChVector<> vz(0, 0, 1);

    // Output in order
    outstream << pos.x() << ",";
    outstream << pos.y() << ",";
    outstream << pos.z() << ",";
    outstream << vx.x() << ",";
    outstream << vx.y() << ",";
    outstream << vx.z() << ",";
    outstream << vy.x() << ",";
    outstream << vy.y() << ",";
    outstream << vy.z() << ",";
    outstream << vz.x() << ",";
    outstream << vz.y() << ",";
    outstream << vz.z() << ",";
    outstream << mesh_scaling.x << "," << mesh_scaling.y << "," << mesh_scaling.z;
    outstream << "\n";
}

int main(int argc, char* argv[]) {
    granular::SetDataPath(std::string(PROJECTS_DATA_DIR) + "granular/");
    sim_param_holder params;

    // Some of the default values might be overwritten by user via command line
    if (argc != 5 || ParseJSON(argv[1], params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }

    float iteration_step = params.step_size;
    const float Bx = params.box_X;
    const float By = params.box_Y;
    const float Bz = params.box_Z;

    // Overwrite parameters from the command line
    params.sphere_radius = std::stof(argv[3]);
    params.sphere_density = std::stof(argv[4]);
    std::cout << "sphere_radius " << params.sphere_radius << std::endl;
    std::cout << "sphere_density " << params.sphere_density << std::endl;

    ChGranularChronoTriMeshAPI apiSMC_TriMesh(params.sphere_radius, params.sphere_density, make_float3(Bx, By, Bz));

    ChSystemGranularSMC_trimesh& gran_sys = apiSMC_TriMesh.getGranSystemSMC_TriMesh();

    gran_sys.setVerbose(params.verbose);

    gran_sys.set_K_n_SPH2SPH(params.normalStiffS2S);
    gran_sys.set_K_n_SPH2WALL(params.normalStiffS2W);
    gran_sys.set_K_n_SPH2MESH(params.normalStiffS2M);

    gran_sys.set_K_t_SPH2SPH(params.tangentStiffS2S);
    gran_sys.set_K_t_SPH2WALL(params.tangentStiffS2W);
    gran_sys.set_K_t_SPH2MESH(params.tangentStiffS2M);

    gran_sys.set_Gamma_n_SPH2SPH(params.normalDampS2S);
    gran_sys.set_Gamma_n_SPH2WALL(params.normalDampS2W);
    gran_sys.set_Gamma_n_SPH2MESH(params.normalDampS2M);

    gran_sys.set_Gamma_t_SPH2SPH(params.tangentDampS2S);
    gran_sys.set_Gamma_t_SPH2WALL(params.tangentDampS2W);
    gran_sys.set_Gamma_t_SPH2MESH(params.tangentDampS2M);

    gran_sys.set_Cohesion_ratio(params.cohesion_ratio);
    gran_sys.set_Adhesion_ratio_S2W(params.adhesion_ratio_s2w);
    gran_sys.set_Adhesion_ratio_S2M(params.adhesion_ratio_s2m);

    gran_sys.set_friction_mode(GRAN_FRICTION_MODE::FRICTIONLESS);  // TODO may want to play with this
    gran_sys.set_static_friction_coeff_SPH2SPH(params.static_friction_coeffS2S);
    gran_sys.set_static_friction_coeff_SPH2WALL(params.static_friction_coeffS2W);
    gran_sys.set_static_friction_coeff_SPH2MESH(params.static_friction_coeffS2M);

    gran_sys.set_rolling_mode(GRAN_ROLLING_MODE::NO_RESISTANCE);  // TODO may want to play with this

    gran_sys.setOutputMode(params.write_mode);
    std::string out_dir(argv[2]);
    filesystem::create_directory(filesystem::path(out_dir));

    gran_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    gran_sys.set_fixed_stepSize(params.step_size);
    gran_sys.set_BD_Fixed(true);
    gran_sys.set_gravitational_acceleration(params.grav_X, params.grav_Y, params.grav_Z);

    // Fill the entire height
    const float fill_bottom = -Bz / 2.f;
    const float fill_height = Bz;

    // Cylinder mesh has interior radius 1 and total radius 1.1
    float cyl_center[3] = {0, 0, 0};
    const float3 scaling = make_float3(Bx / 4.f, Bx / 4.f, Bz);
    std::cout << "Cylinder radius: " << scaling.x << std::endl;

    utils::PDSampler<float> sampler(2.05 * params.sphere_radius);
    std::vector<ChVector<float>> body_points;

    const float fill_radius = scaling.x - 2.f * params.sphere_radius;
    const float fill_top = fill_bottom + fill_height;
    std::cout << "Fill radius " << fill_radius << std::endl;
    std::cout << "Fill bottom " << fill_bottom << std::endl;
    std::cout << "Fill top " << fill_top << std::endl;

    ChVector<float> center(0, 0, fill_bottom);
    center.z() += 2 * params.sphere_radius;
    while (center.z() < fill_top - 2 * params.sphere_radius) {
        auto points = sampler.SampleCylinderZ(center, fill_radius, 0);
        body_points.insert(body_points.end(), points.begin(), points.end());
        center.z() += 2.05 * params.sphere_radius;
    }

    unsigned int n_spheres = body_points.size();
    std::cout << "Adding " << n_spheres << " particles" << std::endl;
    apiSMC_TriMesh.setElemsPositions(body_points);

    std::vector<string> mesh_filenames;
    std::vector<ChMatrix33<float>> mesh_rotscales;
    std::vector<float3> mesh_translations;
    std::vector<float> mesh_masses;
    const float mass = 10;

    std::string mesh_filename = granular::GetDataFile("test_GRAN_commvessels/cylinder_refined.obj");

    mesh_filenames.push_back(mesh_filename);
    mesh_rotscales.push_back(ChMatrix33<float>(ChVector<float>(scaling.x, scaling.y, scaling.z)));
    mesh_translations.push_back(make_float3(0, 0, 0));
    mesh_masses.push_back(mass);

    std::vector<bool> mesh_inflated;
    mesh_inflated.push_back(false);
    std::vector<float> mesh_inflation_radii;
    mesh_inflation_radii.push_back(0);
    apiSMC_TriMesh.load_meshes(mesh_filenames, mesh_rotscales, mesh_translations, mesh_masses, mesh_inflated,
                               mesh_inflation_radii);

    double* meshPosRot = new double[7];
    float* meshVel = new float[6]();

    gran_sys.initialize();
    std::cout << "Writing init..." << std::endl;
    gran_sys.writeFile(out_dir + std::string("/init"));

    const double time_settling = std::sqrt(-2.0 * (params.box_Z) / params.grav_Z);
    const double raising_vel = 1.0;

    const double raising_dist = 10 * 2.0 * 0.2;  // Hard-coded to be the same height as the 0.2 radius run
    const double time_raising = raising_dist / raising_vel;
    const double time_sitting = 10.0;  // TODO no idea how much is enough

    std::cout << "Time settling " << time_settling << std::endl;
    std::cout << "Raising velocity " << raising_vel << std::endl;
    std::cout << "Raising distance " << raising_dist << std::endl;
    std::cout << "Time raising " << time_raising << std::endl;
    std::cout << "Time sitting " << time_sitting << std::endl;

    double mesh_z = 0.0;
    double mesh_vz = raising_vel;

    // Set initial mesh locations for the settling phase
    meshPosRot[0] = 0;
    meshPosRot[1] = 0;
    meshPosRot[2] = mesh_z;
    meshPosRot[3] = 1;
    meshPosRot[4] = 0;
    meshPosRot[5] = 0;
    meshPosRot[6] = 0;

    gran_sys.meshSoup_applyRigidBodyMotion(meshPosRot, meshVel);

    unsigned int currframe = 0;
    double out_fps = 60;
    float frame_step = 1.f / out_fps;  // Duration of a frame
    unsigned int out_steps = frame_step / iteration_step;
    std::cout << "Writing at " << out_fps << " FPS" << std::endl;

    unsigned int step = 0;
    bool settled = false;
    bool raised = false;
    ChVector<> pos_mesh(0, 0, mesh_z);

    std::cout << "Settling..." << std::endl;
    for (float t = 0; t < time_settling + time_raising + time_sitting; t += iteration_step, step++) {
        if (t >= time_settling && t <= time_settling + time_raising) {
            // Raising phase
            if (!settled) {
                std::cout << "Raising..." << std::endl;
                settled = true;
            }

            mesh_z += iteration_step * raising_vel;
            meshPosRot[2] = mesh_z;
            pos_mesh.z() = mesh_z;  // For output
            meshVel[2] = mesh_vz;
            gran_sys.meshSoup_applyRigidBodyMotion(meshPosRot, meshVel);
        } else if (t > time_settling + time_raising) {
            if (!raised) {
                std::cout << "Raised." << std::endl;
                raised = true;
                meshVel[2] = 0;
                gran_sys.meshSoup_applyRigidBodyMotion(meshPosRot, meshVel);
            }
        }

        if (step % out_steps == 0) {
            std::cout << "Rendering frame " << currframe << std::endl;
            char filename[100];
            sprintf(filename, "%s/step%06u", out_dir.c_str(), currframe++);
            gran_sys.writeFile(std::string(filename));
            // gran_sys.write_meshes(std::string(filename));
            std::string mesh_output = std::string(filename) + "_meshframes.csv";

            std::ofstream meshfile(mesh_output);
            std::ostringstream outstream;
            outstream << "mesh_name,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3,sx,sy,sz\n";

            writeMeshFrames(outstream, mesh_filename, pos_mesh, scaling);

            meshfile << outstream.str();
            meshfile.close();
        }

        gran_sys.advance_simulation(iteration_step);
    }

    delete[] meshPosRot;
    delete[] meshVel;
    return 0;
}
