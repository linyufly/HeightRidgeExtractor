// Author: Mingcheng Chen (linyufly@gmail.com)

#include "height_ridge_extractor.h"

#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkIndent.h>

#include <cstdio>
#include <cstdlib>

#include <iostream>

const char *kDataFile = "data/gyre_half.vtk";
const int kNumberOfSamples = 10;

void get_gradient_and_hessian_test() {
  printf("get_gradient_and_hessian_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();

  reader->SetFileName(kDataFile);
  reader->Update();

  vtkStructuredPoints *gradient = NULL, *hessian = NULL;
  HeightRidgeExtractor::get_gradient_and_hessian(reader->GetOutput(),
                                                 &gradient,
                                                 &hessian);

  gradient->PrintSelf(std::cout, vtkIndent(0));
  hessian->PrintSelf(std::cout, vtkIndent(0));

  double tensor[9];
  for (int i = 0; i < kNumberOfSamples; i++) {
    hessian->GetPointData()->GetScalars()->GetTuple(i, tensor);
    for (int row = 0; row < 3; row++) {
      for (int col = 0; col < 3; col++) {
        if (col) printf(" ");
        printf("%lf", tensor[row * 3 + col]);
      }
      printf("\n");
    }
    printf("\n");
  }

  printf("} get_gradient_and_hessian_test\n\n");
}

int main() {
  get_gradient_and_hessian_test();

  return 0;
}
