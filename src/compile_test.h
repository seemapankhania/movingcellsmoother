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

  // Define my custom cell MyCell, which extends Cell by adding extra data members: cell_colour
BDM_SIM_OBJECT(MyCell, bdm::Cell) { // our object extends the Cell object
  BDM_SIM_OBJECT_HEADER(MyCellExt, 1, cell_color_, prev_movement_); // create the header with our new data member

  public:
    MyCellExt() {}
    MyCellExt(const std::array<double, 3>& position) : Base(position) {} // our creator
    // getter and setter for our new data member
    void SetCellColor(int cellcolor) { cell_color_[kIdx] = cellcolor; }
    int GetCellColor() { return cell_color_[kIdx]; }
    int* GetCellColorPtr() { return cell_color_.data(); }

    void SavePositionUpdate(std::array<double, 3> prev_movement) {prev_movement_[kIdx] = prev_movement; }
    std::array<double,3> GetPositionUpdate() { return prev_movement_[kIdx]; }

  private:
  // private data can only be accessed by public function and not directly
    vec<int> cell_color_; // declare our new data member and define its type
    vec<std::array<double,3>> prev_movement_;
};


using namespace std;

struct GrowthModule : public BaseBiologyModule {
  GrowthModule() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {
    auto* random = TSimulation::GetActive()->GetRandom();


    //if (cell->GetDiameter() < 8) {
          //cell->ChangeVolume(400);

          array<double, 3> cell_movements{random->Uniform(-20, 20), random->Uniform(-20, 20), random->Uniform(-20, 20)}; // create an array of 3 random numbers between -2 and 2


          //cell->UpdatePosition(cell_movements); // update the cell mass location, ie move the cell

          array<double, 3> last_movement = cell->GetPositionUpdate();
          array<double, 3> curr_movement = {0.5*last_movement[0]+0.02*cell_movements[0],0.5*last_movement[1]+0.02*cell_movements[1],0.5*last_movement[2]+0.02*cell_movements[2]};
          cell->SavePositionUpdate(curr_movement);
          cell->UpdatePosition(curr_movement);


          cell->SavePositionUpdate(cell_movements);
          cell->SetPosition(cell->GetPosition()); // set the cell position
          cell->SetCellcolour(cell->GetCellcolour());

      } // end Run

  //ClassDefNV(GrowthModule, 1);
};

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<GrowthModule>;  // add GrowthModule
  using AtomicTypes = VariadicTypedef<MyCell>;   // use MyCell object
};

inline int Simulate(int argc, const char** argv) {
  Simulation<> simulation(argc, argv);
  auto* rm = simulation.GetResourceManager();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  #pragma omp parallel
   simulation.GetRandom()->SetSeed(omp_get_thread_num()); //these 2 lines initialise the random number generator so that the cells have different thread values
  size_t nb_of_cells = 10;  // number of cells in the simulation
  double x_coord, y_coord, z_coord;

  param->bound_space_ = true;
  param->min_bound_ = 0;
  param->max_bound_ = 100;  // cube of 100*100*100
  // create a structure to contain cells
  auto* cells = rm->template Get<MyCell>();
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
    MyCell cell({x_coord, y_coord, z_coord});
    // set cell parameters
    cell.SetDiameter(7);
    cell.AddBiologyModule(GrowthModule());
    cell.SetCellColor(4);
    cell.SavePositionUpdate({0.0,0.0,0.0});

    //cell.SetCellColor(static_cast<int>((y_coord / param->max_bound_ * 6)));
    // will vary from 0 to 5. so 6 different layers depending on y_coord
    // cell.SetCellcolour((int)(y_coord / param->max_bound_ * 6));

    cells->push_back(cell);  // put the created cell in our cells structure
  }
  cells->Commit();
  // Run simulation for one timestep
  simulation.GetScheduler()->Simulate(300);

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // COMPILE_TEST_H_
