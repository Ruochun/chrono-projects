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
// Authors: Alessandro Tasora
// =============================================================================
//
// FEA for 3D beams and constrains using the Matlab engine
//
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkLock.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkRackpinion.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMesh.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include "chrono_matlab/ChMatlabEngine.h"
#include "chrono_matlab/ChSolverMatlab.h"

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

int main(int argc, char* argv[]) {
    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    // Create a Chrono::Engine physical system
    ChSystemNSC sys;

    double scales = 100;


    double thickZ = scales * 0.00015;
    double hbarW = scales * 0.00070;
    double hbarL1 = scales * 0.00381;
    double hbarL2 = scales * 0.00387;
    double hbarL3 = scales * 0.00381;
    double vbarL = scales * 0.01137;
    double vbarW = scales * 0.00006;
    double Rpinion = scales * 0.00040;
    double OffPin = scales * 0.00050;
    double Rbalance = scales * 0.00500;
    double Wbalance = scales * 0.00015;
    bool simple_rack = false;

    ChVector<> vAh(-hbarL1 - hbarL2 * 0.5, vbarL, 0);
    ChVector<> vBh(-hbarL2 * 0.5, vbarL, 0);
    ChVector<> vCh(hbarL2 * 0.5, vbarL, 0);
    ChVector<> vDh(hbarL1 + hbarL2 * 0.5, vbarL, 0);
    ChVector<> vAl(-hbarL1 - hbarL2 * 0.5, 0, 0);
    ChVector<> vBl(-hbarL2 * 0.5, 0, 0);
    ChVector<> vCl(hbarL2 * 0.5, 0, 0);
    ChVector<> vDl(hbarL1 + hbarL2 * 0.5, 0, 0);
    ChVector<> vP(0, -Rpinion - hbarW * 0.5, 0);

    // Create a truss:
    auto body_truss = chrono_types::make_shared<ChBody>();

    body_truss->SetBodyFixed(true);

    sys.AddBody(body_truss);

    /*
    // Attach a 'box' shape asset for visualization.
    auto mboxtruss = chrono_types::make_shared<ChBoxShape>();
    mboxtruss->GetBoxGeometry().Pos  = ChVector<>(-0.01, 0,0);
    mboxtruss->GetBoxGeometry().SetLengths( ChVector<>(0.02, 0.2, 0.1) );
    body_truss->AddAsset(mboxtruss);
    */

    // Create a FEM mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create the horizontal beams

    auto msectionH = chrono_types::make_shared<ChBeamSectionAdvanced>();

    msectionH->SetDensity(7000);  //***TEST*** must be 7k
    msectionH->SetYoungModulus(200.0e9);
    msectionH->SetGwithPoissonRatio(0.32);
    msectionH->SetAsRectangularSection(hbarW, thickZ);
    msectionH->SetBeamRaleyghDamping(0.00);

    ChBuilderBeamEuler builder;

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionH,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      2,                     // the number of ChElementBeamEuler to create
                      vAh,                   // the 'Ah' point in space (beginning of beam)
                      vBh,                   // the 'Bh' point in space (end of beam)
                      ChVector<>(0, 1, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Ah = builder.GetLastBeamNodes().front();
    std::shared_ptr<ChNodeFEAxyzrot> node_Bh = builder.GetLastBeamNodes().back();

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionH,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      2,                     // the number of ChElementBeamEuler to create
                      node_Bh,               // the 'Bh' point in space (beginning of beam)
                      vCh,                   // the 'Ch' point in space (end of beam)
                      ChVector<>(0, 1, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Ch = builder.GetLastBeamNodes().back();

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionH,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      2,                     // the number of ChElementBeamEuler to create
                      node_Ch,               // the 'Ch' point in space (beginning of beam)
                      vDh,                   // the 'Dh' point in space (end of beam)
                      ChVector<>(0, 1, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Dh = builder.GetLastBeamNodes().back();

    // Create the vertical flexible beams

    auto msectionV = chrono_types::make_shared<ChBeamSectionAdvanced>();

    msectionV->SetDensity(7000);  //***TEST*** must be 7k
    msectionV->SetYoungModulus(200.0e9);
    msectionV->SetGwithPoissonRatio(0.32);
    msectionV->SetAsRectangularSection(vbarW, thickZ);
    msectionV->SetBeamRaleyghDamping(0.00);

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionV,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      6,                     // the number of ChElementBeamEuler to create
                      node_Ah,               // the 'Ah' point in space (beginning of beam)
                      vAl,                   // the 'Al' point in space (end of beam)
                      ChVector<>(1, 0, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Al = builder.GetLastBeamNodes().back();

    node_Al->SetFixed(true);

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionV,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      6,                     // the number of ChElementBeamEuler to create
                      node_Dh,               // the 'Dh' point in space (beginning of beam)
                      vDl,                   // the 'Dl' point in space (end of beam)
                      ChVector<>(1, 0, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Dl = builder.GetLastBeamNodes().back();

    node_Dl->SetFixed(true);

    // Create the inner vertical flexible beams

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionV,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      6,                     // the number of ChElementBeamEuler to create
                      node_Bh,               // the 'Bh' point in space (beginning of beam)
                      vBl,                   // the 'Bl' point in space (end of beam)
                      ChVector<>(1, 0, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Bl = builder.GetLastBeamNodes().back();

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionV,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      6,                     // the number of ChElementBeamEuler to create
                      node_Ch,               // the 'Dh' point in space (beginning of beam)
                      vCl,                   // the 'Dl' point in space (end of beam)
                      ChVector<>(1, 0, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Cl = builder.GetLastBeamNodes().back();

    // Create the rack

    if (simple_rack) {
        builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                          msectionH,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                          1,                     // the number of ChElementBeamEuler to create
                          node_Bl,               // the 'Cl' point in space (beginning of beam)
                          node_Cl,               // the 'Dl' point in space (end of beam)
                          ChVector<>(0, 1, 0));  // the 'Y' up direction of the section for the beam
    }

    //
    // Final touches to mesh..
    //

    // Remember to add the mesh to the system!
    sys.Add(my_mesh);

    // ==Asset== attach a visualization of the FEM mesh.
    // This will automatically update a triangle mesh (a ChTriangleMeshShape
    // asset that is internally managed) by setting  proper
    // coordinates and vertex colours as in the FEM elements.
    // Such triangle mesh can be rendered by Irrlicht or POVray or whatever
    // postprocessor that can handle a coloured ChTriangleMeshShape).
    // Do not forget AddAsset() at the end!

    auto mvisualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    mvisualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_SPEED_NORM);
    mvisualizebeamA->SetColorscaleMinMax(-30, 30);
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    my_mesh->AddVisualShapeFEA(mvisualizebeamA);

    auto mvisualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    mvisualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_CSYS);
    mvisualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    mvisualizebeamC->SetSymbolsThickness(0.001);
    mvisualizebeamC->SetSymbolsScale(0.01);
    mvisualizebeamC->SetZbufferHide(false);
    my_mesh->AddVisualShapeFEA(mvisualizebeamC);

    //
    // The balance and the rigid rach
    //

    if (!simple_rack) {
        auto rack = chrono_types::make_shared<ChBodyEasyBox>(hbarL2, hbarW, thickZ, 7000, true, false);
        rack->SetPos(0.5 * (vBl + vCl));
        sys.Add(rack);

        auto constr_B = chrono_types::make_shared<ChLinkMateGeneric>();
        constr_B->Initialize(node_Bl, rack, false, node_Bl->Frame(), node_Bl->Frame());
        sys.Add(constr_B);

        auto constr_C = chrono_types::make_shared<ChLinkMateGeneric>();
        constr_C->Initialize(node_Cl, rack, false, node_Cl->Frame(), node_Cl->Frame());
        sys.Add(constr_C);

        auto balance = chrono_types::make_shared<ChBodyEasyCylinder>(Rbalance, Wbalance, 7000, true, false);
        balance->SetPos(vP + ChVector<>(0, 0, -OffPin));
        balance->SetRot(Q_from_AngAxis(CH_C_PI_2, VECT_X));
        for (int i = 0; i < 6; ++i) {
            double phi = CH_C_2PI * (i / 6.0);
            auto vshape = chrono_types::make_shared<ChCylinderShape>();
            vshape->GetCylinderGeometry().p1 =
                ChVector<>(sin(phi) * Rbalance * 0.8, Wbalance, cos(phi) * Rbalance * 0.8);
            vshape->GetCylinderGeometry().p2 = vshape->GetCylinderGeometry().p1 + ChVector<>(0, 2 * Wbalance, 0);
            vshape->GetCylinderGeometry().rad = Rbalance * 0.1;
            balance->AddVisualShape(vshape);
        }
        auto vshaft = chrono_types::make_shared<ChCylinderShape>();
        vshaft->GetCylinderGeometry().p1 = vP + ChVector<>(0, -OffPin * 10, 0);
        vshaft->GetCylinderGeometry().p2 = vP + ChVector<>(0, OffPin * 10, 0);
        vshaft->GetCylinderGeometry().rad = Rpinion;
        vshaft->SetColor(ChColor(0.5f, 0.9f, 0.9f));
        balance->AddVisualShape(vshaft);

        sys.Add(balance);

        auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
        std::shared_ptr<ChBody> mbalance = balance;
        revolute->Initialize(mbalance, body_truss, ChCoordsys<>(vP + ChVector<>(0, 0, -0.01)));

        sys.Add(revolute);

        auto constr_rack = chrono_types::make_shared<ChLinkRackpinion>();
        constr_rack->Initialize(balance, rack, false, ChFrame<>(), ChFrame<>());

        ChFrameMoving<> f_pin_abs(vP);
        ChFrameMoving<> f_rack_abs(vP + ChVector<>(0, 0.1, 0));
        ChFrameMoving<> f_pin;
        ChFrameMoving<> f_rack;
        balance->TransformParentToLocal(f_pin_abs, f_pin);
        rack->TransformParentToLocal(f_rack_abs, f_rack);
        constr_rack->SetPinionRadius(Rpinion);
        constr_rack->SetPinionFrame(f_pin);
        constr_rack->SetRackFrame(f_rack);

        sys.Add(constr_rack);

        balance->SetWvel_par(ChVector<>(0, 0, 1.5));
    }
    
    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Beams and constraints");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(0, 0.01*scales, 0.01*scales));
    vis->AddTypicalLights();
    vis->AttachSystem(&sys);

    //
    // THE SOFT-REAL-TIME CYCLE
    //

    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    solver->EnableDiagonalPreconditioner(false);
    solver->EnableWarmStart(true);
    solver->SetMaxIterations(400);
    solver->SetVerbose(true);
    sys.SetSolver(solver);
    sys.SetSolverForceTolerance(1e-12);

    //***TEST***
    ChMatlabEngine matlab_engine;
    auto matlab_solver = chrono_types::make_shared<ChSolverMatlab>(matlab_engine);
    sys.SetSolver(matlab_solver);

    sys.Set_G_acc(ChVector<>(0, 0, 0));

    GetLog() << "STATIC linear solve ----\n";
    node_Cl->SetForce(ChVector<>(50, 0, 0));
    // application.GetSystem()->DoStaticLinear();
    node_Cl->SetForce(ChVector<>(0, 0, 0));

    if (simple_rack) {
        node_Cl->SetForce(ChVector<>(50, 0, 0));
        sys.DoStaticNonlinear(12);
        node_Cl->SetForce(ChVector<>(0, 0, 0));
    }

    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        tools::drawGrid(vis.get(), 0.2, 0.2, 10, 10, ChCoordsys<>(VNULL, CH_C_PI_2, VECT_Z), ChColor(0.4f, 0.5f, 0.5f),
                        true);
        vis->EndScene();
        sys.DoStepDynamics(0.01);
    }

    return 0;
}
