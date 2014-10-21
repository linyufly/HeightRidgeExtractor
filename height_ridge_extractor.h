// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef HEIGHT_RIDGE_EXTRACTOR_H_
#define HEIGHT_RIDGE_EXTRACTOR_H_

class vtkStructuredPoints;
class vtkPolyData;

class HeightRidgeExtractor {
 public:
   vtkPolyData *extract_ridges(const vtkStructuredPoints &scalar_field);
}

#endif  // HEIGHT_RIDGE_EXTRACTOR_H_
