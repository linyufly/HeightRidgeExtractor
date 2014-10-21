// Author: Mingcheng Chen (linyufly@gmail.com)

#include "height_ridge_extractor.h"

#include <vtkStructuredPoints.h>
#include <vtkPolyData.h>
#include <vtkGradientFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

void HeightRidgeExtractor::get_gradient_and_hessian(
    vtkStructuredPoints *scalar_field,
    vtkStructuredPoints **gradient,
    vtkStructuredPoints **hessian) {
  // Calculate the gradient
  vtkSmartPointer<vtkGradientFilter> grad_filter =
      vtkSmartPointer<vtkGradientFilter>::New();
  grad_filter->SetInputData(scalar_field);
  grad_filter->Update();

  *gradient = vtkStructuredPoints::New();
  (*gradient)->ShallowCopy(grad_filter->GetOutput());

  (*gradient)->GetPointData()->SetActiveScalars("Gradients");
  (*gradient)->GetPointData()->GetScalars()->SetName("gradient");

  // Calculate the hessian
  grad_filter = vtkSmartPointer<vtkGradientFilter>::New();
  grad_filter->SetInputData(*gradient);
  grad_filter->Update();

  *hessian = vtkStructuredPoints::New();
  (*hessian)->ShallowCopy(grad_filter->GetOutput());

  (*hessian)->GetPointData()->SetActiveScalars("Gradients");
  (*hessian)->GetPointData()->GetScalars()->SetName("hessian");
}

vtkPolyData *HeightRidgeExtractor::extract_ridges(
    vtkStructuredPoints *scalar_field) {
  // Not implemented
}
