#include "helpers/benchmark_result.h"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/ipc/feather.h>

#include <iostream>
#include <string>

#include "helpers/arrow_utils.h"

namespace bench {

namespace {

arrow::Status SaveToFeatherFile(const BenchmarkResultArray& res_arr,
                                const std::string& filename) {
  std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("solve_time", arrow::float64()),
      arrow::field("prim_dual_gap", arrow::float64()),
      arrow::field("prim_infeas_err", arrow::float64()),
      arrow::field("dual_infeas_err", arrow::float64()),
      arrow::field("iter", arrow::float64()),
      arrow::field("optimal", arrow::boolean())};

  if (res_arr.store_polytope_size) {
    fields.push_back(arrow::field("polytope_size", arrow::int32()));
  }

  std::shared_ptr<arrow::Schema> schema = arrow::schema(fields);

  std::shared_ptr<arrow::Array> solve_time_array;
  ARROW_ASSIGN_OR_RAISE(solve_time_array, MakeDoubleArray(res_arr.solve_times));

  std::shared_ptr<arrow::Array> prim_dual_gap_array;
  ARROW_ASSIGN_OR_RAISE(prim_dual_gap_array,
                        MakeDoubleArray(res_arr.prim_dual_gaps));

  std::shared_ptr<arrow::Array> prim_infeas_err_array;
  ARROW_ASSIGN_OR_RAISE(prim_infeas_err_array,
                        MakeDoubleArray(res_arr.prim_infeas_errs));

  std::shared_ptr<arrow::Array> dual_infeas_err_array;
  ARROW_ASSIGN_OR_RAISE(dual_infeas_err_array,
                        MakeDoubleArray(res_arr.dual_infeas_errs));

  std::shared_ptr<arrow::Array> iter_array;
  ARROW_ASSIGN_OR_RAISE(iter_array, MakeDoubleArray(res_arr.iters));

  std::shared_ptr<arrow::Array> optimal_flag_array;
  ARROW_ASSIGN_OR_RAISE(optimal_flag_array,
                        MakeBooleanArray(res_arr.optimal_flags));

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      solve_time_array,      prim_dual_gap_array, prim_infeas_err_array,
      dual_infeas_err_array, iter_array,          optimal_flag_array};

  if (res_arr.store_polytope_size) {
    std::shared_ptr<arrow::Array> polytope_size_array;
    ARROW_ASSIGN_OR_RAISE(polytope_size_array,
                          MakeInt32Array(res_arr.polytope_sizes));
    arrays.push_back(polytope_size_array);
  }

  std::shared_ptr<arrow::Table> table = arrow::Table::Make(schema, arrays);

  std::shared_ptr<arrow::io::FileOutputStream> file;
  ARROW_ASSIGN_OR_RAISE(file, arrow::io::FileOutputStream::Open(filename));

  ARROW_RETURN_NOT_OK(arrow::ipc::feather::WriteTable(*table, file.get()));

  ARROW_RETURN_NOT_OK(file->Close());
  return arrow::Status::OK();
}

}  // namespace

bool BenchmarkResultArray::SaveToFile(const std::string& filename) const {
  arrow::Status status = SaveToFeatherFile(*this, filename);
  if (!status.ok()) {
    std::cerr << "Failed to write Feather file: " << status.ToString()
              << std::endl;
    return false;
  }

  std::cout << "Benchmark data written to: " << filename << std::endl;
  return true;
}

}  // namespace bench
