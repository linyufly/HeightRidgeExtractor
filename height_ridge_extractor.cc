// Author: Mingcheng Chen (linyufly@gmail.com)

#include "height_ridge_extractor.h"
#include "util.h"
#include "math.h"

#include <vtkStructuredPoints.h>
#include <vtkPolyData.h>
#include <vtkGradientFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <cstring>

namespace {

const int kVertexList[8][3] = {
    {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
};

const int kEdgeList[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
};

double dot_product_3d(double *a, double *b) {
  double result = 0.0;
  for (int i = 0; i < 3; i++) {
    result += a[i] * b[i];
  }
  return result;
}

}

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
  int dimensions[3];
  double spacing[3], origin[3];
  scalar_field->GetDimensions(dimensions);
  scalar_field->GetSpacing(spacing);
  scalar_field->GetOrigin(origin);

  // Calculate the gradient and hessian
  vtkStructuredPoints *gradient_field = NULL, *hessian_field = NULL;
  get_gradient_and_hessian(scalar_field, &gradient_field, &hessian_field);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int *mark_data = new int[nx * ny * nz * 3];
  memset(mark_data, 255, sizeof(int) * nx * ny * nz * 3);

  int **mark_xyz = new int *[nx * ny * nz];
  for (int i = 0; i < nx * ny * nz; i++) {
    mark_xyz[i] = mark_data + i * 3;
  }

  int ***mark_xy = new int **[nx * ny];
  for (int i = 0; i < nx * ny; i++) {
    mark_xy[i] = mark_xyz + i * nz;
  }

  // mark_x[i][j][k][d] indicates a unique edge.
  int ****mark_x = new int ***[nx];
  for (int i = 0; i < nx; i++) {
    mark_x[i] = mark_xy + i * ny;
  }

  int num_vertices = 0;

  for (int x = 0; x + 1 < nx; x++) {
    for (int y = 0; y + 1 < ny; y++) {
      for (int z = 0; z + 1 < nz; z++) {
        double dot_prod[3][3][3], e3[3][3][3][3], grad[3][3][3][3];

        // Collect e3 and grad
        for (int dx = 0; dx < 2; dx++) {
          for (int dy = 0; dy < 2; dy++) {
            for (int dz = 0; dz < 2; dz++) {
              int curr_x = x + dx;
              int curr_y = y + dy;
              int curr_z = z + dz;

              double **hessian = create_matrix<double>(3, 3);
              int point_id = (curr_z * ny + curr_y) * nx + curr_x;

              double tensor[9];
              hessian_field->GetPointData()->GetScalars()
                                           ->GetTuple(point_id, tensor);
              for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                  hessian[i][j] = tensor[i * 3 + j];
                }
              }

              double *eigen_values = new double[3];
              double **eigen_vectors = create_matrix<double>(3, 3);
              vtkMath::Jacobi(hessian, eigen_values, eigen_vectors);

              // e3 indicates the eigen vector with the smallest eigen value.
              for (int i = 0; i < 3; i++) {
                e3[dx][dy][dz][i] = eigen_vectors[i][2];
              }

              delete [] eigen_values;
              delete_matrix(eigen_vectors);
              delete_matrix(hessian);

              gradient_field->GetPointData()->GetScalars()
                                            ->GetTuple(point_id, tensor);
              for (int i = 0; i < 3; i++) {
                grad[dx][dy][dz] = tensor[i];
              }
            }
          }
        }

        // Re-orientate e3
        double **vectors = create_matrix<double>(8, 3);
        int num_vectors = 0;
        for (int dx = 0; dx < 2; dx++) {
          for (int dy = 0; dy < 2; dy++) {
            for (int dz = 0; dz < 2; dz++) {
              for (int i = 0; i < 3; i++) {
                vectors[num_vectors][i] = e3[dx][dy][dz][i];
              }
              num_vectors++;
            }
          }
        }

        double *pivot = principal_component(vectors, 8, 3);

        delete_matrix(vectors);

        for (int dx = 0; dx < 2; dx++) {
          for (int dy = 0; dy < 2; dy++) {
            for (int dz = 0; dz < 2; dz++) {
              if (dot_product_3d(pivot, e3[dx][dy][dz]) < 0.0) {
                for (int i = 0; i < 3; i++) {
                  e3[dx][dy][dz][i] *= -1.0;
                }
              }

              dot_prod[dx][dy][dz] = dot_product_3d(grad[dx][dy][dz],
                                                    e3[dx][dy][dz]);
            }
          }
        }

        delete [] pivot;
      }
    }
  }

  gradient_field->Delete();
  hessian_field->Delete();

  return NULL;
}
