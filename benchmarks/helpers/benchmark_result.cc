#include "helpers/benchmark_result.h"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/ipc/feather.h>

#include <iostream>
#include <string>
#include <vector>

namespace internal {

namespace {

arrow::Result<std::shared_ptr<arrow::Array>> MakeDoubleArray(
    const std::vector<double>& data) {
  arrow::DoubleBuilder builder;
  ARROW_RETURN_NOT_OK(builder.AppendValues(data.data(), data.size()));
  std::shared_ptr<arrow::Array> array;
  ARROW_RETURN_NOT_OK(builder.Finish(&array));
  return array;
}

// Assumes int maps to arrow::int32().
arrow::Result<std::shared_ptr<arrow::Array>> MakeInt32Array(
    const std::vector<int>& data) {
  arrow::Int32Builder builder;
  ARROW_RETURN_NOT_OK(builder.AppendValues(data.data(), data.size()));
  std::shared_ptr<arrow::Array> array;
  ARROW_RETURN_NOT_OK(builder.Finish(&array));
  return array;
}

// Specialization for std::vector<bool> because of its packed nature.
arrow::Result<std::shared_ptr<arrow::Array>> MakeBooleanArray(
    const std::vector<bool>& data) {
  arrow::BooleanBuilder builder;
  ARROW_RETURN_NOT_OK(builder.Reserve(data.size()));
  for (bool val : data) {
    ARROW_RETURN_NOT_OK(builder.Append(val));
  }
  std::shared_ptr<arrow::Array> array;
  ARROW_RETURN_NOT_OK(builder.Finish(&array));
  return array;
}

arrow::Status SaveToFeatherFile(const BenchmarkResultArray& res_arr,
                                const std::string& filename) {
  std::shared_ptr<arrow::Schema> schema =
      arrow::schema({arrow::field("solve_time", arrow::float64()),
                     arrow::field("prim_dual_gap", arrow::float64()),
                     arrow::field("prim_feas_err", arrow::float64()),
                     arrow::field("dual_feas_err", arrow::float64()),
                     arrow::field("iter", arrow::int32()),
                     arrow::field("optimal", arrow::boolean())});

  std::shared_ptr<arrow::Array> solve_time_array;
  ARROW_ASSIGN_OR_RAISE(solve_time_array, MakeDoubleArray(res_arr.solve_times));

  std::shared_ptr<arrow::Array> prim_dual_gap_array;
  ARROW_ASSIGN_OR_RAISE(prim_dual_gap_array,
                        MakeDoubleArray(res_arr.prim_dual_gaps));

  std::shared_ptr<arrow::Array> prim_feas_err_array;
  ARROW_ASSIGN_OR_RAISE(prim_feas_err_array,
                        MakeDoubleArray(res_arr.prim_feas_errs));

  std::shared_ptr<arrow::Array> dual_feas_err_array;
  ARROW_ASSIGN_OR_RAISE(dual_feas_err_array,
                        MakeDoubleArray(res_arr.dual_feas_errs));

  std::shared_ptr<arrow::Array> iter_array;
  ARROW_ASSIGN_OR_RAISE(iter_array, MakeInt32Array(res_arr.iters));

  std::shared_ptr<arrow::Array> optimal_flag_array;
  ARROW_ASSIGN_OR_RAISE(optimal_flag_array,
                        MakeBooleanArray(res_arr.optimal_flags));

  std::shared_ptr<arrow::Table> table{
      arrow::Table::Make(schema, {
                                     solve_time_array,
                                     prim_dual_gap_array,
                                     prim_feas_err_array,
                                     dual_feas_err_array,
                                     iter_array,
                                     optimal_flag_array,
                                 })};

  std::shared_ptr<arrow::io::FileOutputStream> file;
  ARROW_ASSIGN_OR_RAISE(file, arrow::io::FileOutputStream::Open(filename));

  ARROW_RETURN_NOT_OK(arrow::ipc::feather::WriteTable(*table, file.get()));

  ARROW_RETURN_NOT_OK(file->Close());
  return arrow::Status::OK();
}

}  // namespace

bool BenchmarkResultArray::SaveToFile(const std::string& filename) {
  arrow::Status status = SaveToFeatherFile(*this, filename);
  if (!status.ok()) {
    std::cerr << "Failed to write Feather file: " << status.ToString()
              << std::endl;
    return false;
  }

  std::cout << "Benchmark data written to: " << filename << std::endl;
  return true;
}

}  // namespace internal
