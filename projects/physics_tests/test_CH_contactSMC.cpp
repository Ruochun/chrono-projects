//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

#include "chrono/ChConfig.h"
#include "chrono/physics/ChContactContainerSMC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#ifdef CHRONO_PARDISO_MKL
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#endif

#include <irrlicht.h>

using namespace chrono;
using namespace chrono::irrlicht;

// Custom contact container -- get access to the contact lists in the base class.
class MyContactContainer : public ChContactContainerSMC {
  public:
    MyContactContainer() {}

    // Traverse the list contactlist_6_6
    void ScanContacts(std::shared_ptr<ChBody> ball) {
        auto iter = contactlist_6_6.begin();
        while (iter != contactlist_6_6.end()) {
            ChContactable* objA = (*iter)->GetObjA();
            ChContactable* objB = (*iter)->GetObjB();
            ChVector<> pA = (*iter)->GetContactP1();
            ChVector<> pB = (*iter)->GetContactP2();

            GetLog().SetNumFormat("%12.3e");

            int iball;
            if (objA == ball.get()) {
                iball = 0;
                GetLog() << "pA = " << pA << "\n";
            } else if (objB == ball.get()) {
                iball = 1;
                GetLog() << "pB = " << pB << "\n";
            }

            const ChKblockGeneric* KRM = (*iter)->GetJacobianKRM();
            const ChMatrixDynamic<double>* K = (*iter)->GetJacobianK();
            const ChMatrixDynamic<double>* R = (*iter)->GetJacobianR();

            if (KRM) {
                GetLog() << "K = " << *K << "\n";
                GetLog() << "R = " << *R << "\n";

                ChMatrixNM<double, 6, 6> Kball;
                ChMatrixNM<double, 6, 6> Rball;
                Kball = K->block(iball * 6, iball * 6, 6, 6);
                Rball = R->block(iball * 6, iball * 6, 6, 6);

                GetLog() << "Kball = " << Kball << "\n";
                GetLog() << "Rball = " << Rball << "\n";
            }

            ++iter;
        }
    }
};

// ====================================================================================

int main(int argc, char* argv[]) {
    // ---------------------------------
    // Set path to Chrono data directory
    // ---------------------------------
    SetChronoDataPath(CHRONO_DATA_DIR);

    // ---------------------
    // Simulation parameters
    // ---------------------

    double gravity = -9.81;   // gravitational acceleration
    double time_step = 1e-4;  // integration step size

    ChSolver::Type solver_type = ChSolver::Type::PARDISO_MKL;

    bool stiff_contact = true;

    // ---------------------------
    // Contact material properties
    // ---------------------------

    bool use_mat_properties = false;
    ChSystemSMC::ContactForceModel force_model = ChSystemSMC::Hooke;
    ChSystemSMC::TangentialDisplacementModel tdispl_model = ChSystemSMC::None;

    float young_modulus = 2e9f;
    float friction = 0.4f;
    float restitution = 0.1f;
    float adhesion = 0.0f;

    float kn = 1e8;
    float gn = 0;
    float kt = 0;
    float gt = 0;

    // -------------------------------
    // Parameters for the falling ball
    // -------------------------------

    int ballId = 100;
    double radius = 1;
    double mass = 1000;
    ChVector<> pos(0, 2, 0);
    ChQuaternion<> rot(1, 0, 0, 0);
    ChVector<> init_vel(0, 0, 0);
    ChVector<> init_omg(0, 0, 0);

    // ---------------------------------
    // Parameters for the containing bin
    // ---------------------------------

    int binId = 200;
    double width = 2;
    double length = 2;
    double thickness = 0.1;

    // -----------------
    // Create the system
    // -----------------

    ChSystemSMC system(use_mat_properties);

    // Set the SMC contact force model
    system.SetContactForceModel(force_model);

    // Set tangential displacement model
    system.SetTangentialDisplacementModel(tdispl_model);

    // Set contact forces as stiff (to force Jacobian computation) or non-stiff
    system.SetStiffContact(stiff_contact);

    system.Set_G_acc(ChVector<>(0, gravity, 0));

    // Create a material (will be used by both objects)
    auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    material->SetYoungModulus(young_modulus);
    material->SetRestitution(restitution);
    material->SetFriction(friction);
    material->SetAdhesion(adhesion);
    material->SetKn(kn);
    material->SetGn(gn);
    material->SetKt(kt);
    material->SetGt(gt);

    // Create the falling ball
    auto ball = chrono_types::make_shared<ChBody>();

    ball->SetIdentifier(ballId);
    ball->SetMass(mass);
    ball->SetInertiaXX(0.4 * mass * radius * radius * ChVector<>(1, 1, 1));
    ball->SetPos(pos);
    ball->SetRot(rot);
    ball->SetPos_dt(init_vel);
    ball->SetWvel_par(init_omg);
    ball->SetCollide(true);
    ball->SetBodyFixed(false);

    ball->GetCollisionModel()->ClearModel();
    ball->GetCollisionModel()->AddSphere(material, radius);
    ball->GetCollisionModel()->BuildModel();

    auto sphere = chrono_types::make_shared<ChSphereShape>();
    sphere->GetSphereGeometry().rad = radius;
    sphere->SetTexture(GetChronoDataFile("textures/bluewhite.png"));
    ball->AddVisualShape(sphere);

    system.AddBody(ball);

    // Create ground
    auto ground = chrono_types::make_shared<ChBody>();

    ground->SetIdentifier(binId);
    ground->SetMass(1);
    ground->SetPos(ChVector<>(0, 0, 0));
    ground->SetRot(ChQuaternion<>(1, 0, 0, 0));
    ground->SetCollide(true);
    ground->SetBodyFixed(true);

    ground->GetCollisionModel()->ClearModel();
    ground->GetCollisionModel()->AddBox(material, width, thickness, length, ChVector<>(0, -thickness, 0));
    ground->GetCollisionModel()->BuildModel();

    auto box = chrono_types::make_shared<ChBoxShape>();
    box->GetBoxGeometry().Size = ChVector<>(width, thickness, length);
    ground->AddVisualShape(box, ChFrame<>(ChVector<>(0, -thickness, 0)));

    system.AddBody(ground);

    // Create the Irrlicht visualization
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("SMC demo");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(0, 3, -6));
    vis->AddTypicalLights();
    vis->SetSymbolScale(1e-4);
    vis->EnableContactDrawing(ContactsDrawMode::CONTACT_FORCES);
    vis->AttachSystem(&system);

    // ----------------------------
    // Use custom contact container
    // ----------------------------

    auto container = chrono_types::make_shared<MyContactContainer>();
    system.SetContactContainer(container);

    // -------------------
    // Setup linear solver
    // -------------------

    // Note that not all solvers support stiffness matrices (that includes the default SolverSMC).
#ifndef CHRONO_PARDISO_MKL
    if (solver_type == ChSolver::Type::PARDISO_MKL) {
        GetLog() << "PardisoMKL support not enabled.  Solver reset to default.\n";
        solver_type = ChSolver::Type::MINRES;
    }
#endif

    switch (solver_type) {
        case ChSolver::Type::MINRES: {
            GetLog() << "Using MINRES solver.\n";
            auto minres_solver = chrono_types::make_shared<ChSolverMINRES>();
            minres_solver->EnableDiagonalPreconditioner(true);
            minres_solver->SetMaxIterations(100);
            system.SetSolver(minres_solver);
            system.SetSolverForceTolerance(1e-6);
            break;
        }
        case ChSolver::Type::PARDISO_MKL: {
#ifdef CHRONO_PARDISO_MKL
            GetLog() << "Using PardisoMKL solver.\n";
            auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
            mkl_solver->LockSparsityPattern(true);
            system.SetSolver(mkl_solver);
#endif
            break;
        }
        default: {
            GetLog() << "Using DEFAULT solver.\n";
            system.SetSolverMaxIterations(100);
            system.SetSolverForceTolerance(1e-6);
            break;
        }
    }

    // ----------------
    // Setup integrator
    // ----------------

    system.SetTimestepperType(ChTimestepper::Type::HHT);
    auto integrator = std::static_pointer_cast<ChTimestepperHHT>(system.GetTimestepper());
    integrator->SetAlpha(0.0);
    integrator->SetMaxiters(100);
    integrator->SetAbsTolerances(1e-08);
    ////integrator->SetMode(ChTimestepperHHT::POSITION);
    integrator->SetScaling(false);
    ////integrator->SetStepControl(false);
    ////integrator->SetVerbose(true);

    // ---------------
    // Simulation loop
    // ---------------
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        system.DoStepDynamics(time_step);
        if (ball->GetPos().y() <= radius) {
            container->ScanContacts(ball);
            GetLog() << "t = " << system.GetChTime() << "  NR iters. = " << integrator->GetNumIterations() << "\n";
        }
    }

    return 0;
}
