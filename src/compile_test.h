// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef COMPILE_TEST_H_
#define COMPILE_TEST_H_

#include "biodynamo.h"
namespace bdm {

using namespace std;

struct GrowthModule : public BaseBiologyModule {
  GrowthModule() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {
    auto* random = TSimulation::GetActive()->GetRandom();


    if (cell->GetDiameter() < 8) {
          //cell->ChangeVolume(400);

          array<double, 3> cell_movements{random->Uniform(-2, 2), random->Uniform(-2, 2), random->Uniform(-2, 2)}; // create an array of 3 random numbers between -2 and 2
          cell->UpdatePosition(cell_movements); // update the cell mass location, ie move the cell
          // cell->SetPosition(cell->GetPosition()); // set the cell position
        }
      } // end Run

  ClassDefNV(GrowthModule, 1);
};

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<GrowthModule>;  // add GrowthModule
  // using AtomicTypes = VariadicTypedef<MyCell>;   // use MyCell object
};



inline int Simulate(int argc, const char** argv) {
  Simulation<> simulation(argc, argv);
  auto* rm = simulation.GetResourceManager();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  size_t nb_of_cells = 5;  // number of cells in the simulation
  double x_coord, y_coord, z_coord;

  param->bound_space_ = true;
  param->min_bound_ = 0;
  param->max_bound_ = 100;  // cube of 100*100*100

  // create a structure to contain cells
  auto* cells = rm->template Get<Cell>();
  // allocate the correct number of cell in our cells structure before
  // cell creation
  cells->reserve(nb_of_cells);

  for (size_t i = 0; i < nb_of_cells; ++i) {
    // our modelling will be a cell cube of 100*100*100
    // random double between 0 and 100
    x_coord = random->Uniform(param->min_bound_, param->max_bound_);
    y_coord = random->Uniform(param->min_bound_, param->max_bound_);
    z_coord = random->Uniform(param->min_bound_, param->max_bound_);

    // creating the cell at position x, y, z
    Cell cell({x_coord, y_coord, z_coord});
    // set cell parameters
    cell.SetDiameter(7);
    cell.AddBiologyModule(GrowthModule());

    // will vary from 0 to 5. so 6 different layers depending on y_coord
    // cell.SetCellColour((int)(y_coord / param->max_bound_ * 6));

    cells->push_back(cell);  // put the created cell in our cells structure
  }




  // Run simulation for one timestep
  simulation.GetScheduler()->Simulate(50);

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // COMPILE_TEST_H_
