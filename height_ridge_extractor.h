// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef HEIGHT_RIDGE_EXTRACTOR_H_
#define HEIGHT_RIDGE_EXTRACTOR_H_

class vtkStructuredPoints;
class vtkPolyData;

class HeightRidgeExtractor {
 public:
  static void get_gradient_and_hessian(vtkStructuredPoints *scalar_field,
                                       vtkStructuredPoints **gradient,
                                       vtkStructuredPoints **hessian);

  // scalar_field must be 3-dimensional.
  vtkPolyData *extract_ridges(vtkStructuredPoints *scalar_field);
};

#endif  // HEIGHT_RIDGE_EXTRACTOR_H_
