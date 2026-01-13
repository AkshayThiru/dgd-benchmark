#include "helpers/convergence_result.h"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/ipc/feather.h>

#include <iostream>
#include <string>

#include "helpers/arrow_utils.h"

namespace bench {

namespace {

arrow::Status SaveToFeatherFile(const ConvergenceResultArray& res_arr,
                                const std::string& filename) {
  std::shared_ptr<arrow::Schema> schema =
      arrow::schema({arrow::field("min", arrow::float64()),
                     arrow::field("q01", arrow::float64()),
                     arrow::field("q25", arrow::float64()),
                     arrow::field("q50", arrow::float64()),
                     arrow::field("q75", arrow::float64()),
                     arrow::field("q99", arrow::float64()),
                     arrow::field("max", arrow::float64())});

  std::shared_ptr<arrow::Array> min_array;
  ARROW_ASSIGN_OR_RAISE(min_array, MakeDoubleArray(res_arr.min));

  std::shared_ptr<arrow::Array> q01_array;
  ARROW_ASSIGN_OR_RAISE(q01_array, MakeDoubleArray(res_arr.q01));

  std::shared_ptr<arrow::Array> q25_array;
  ARROW_ASSIGN_OR_RAISE(q25_array, MakeDoubleArray(res_arr.q25));

  std::shared_ptr<arrow::Array> q50_array;
  ARROW_ASSIGN_OR_RAISE(q50_array, MakeDoubleArray(res_arr.q50));

  std::shared_ptr<arrow::Array> q75_array;
  ARROW_ASSIGN_OR_RAISE(q75_array, MakeDoubleArray(res_arr.q75));

  std::shared_ptr<arrow::Array> q99_array;
  ARROW_ASSIGN_OR_RAISE(q99_array, MakeDoubleArray(res_arr.q99));

  std::shared_ptr<arrow::Array> max_array;
  ARROW_ASSIGN_OR_RAISE(max_array, MakeDoubleArray(res_arr.max));

  std::shared_ptr<arrow::Table> table{arrow::Table::Make(schema, {
                                                                     min_array,
                                                                     q01_array,
                                                                     q25_array,
                                                                     q50_array,
                                                                     q75_array,
                                                                     q99_array,
                                                                     max_array,
                                                                 })};

  std::shared_ptr<arrow::io::FileOutputStream> file;
  ARROW_ASSIGN_OR_RAISE(file, arrow::io::FileOutputStream::Open(filename));

  ARROW_RETURN_NOT_OK(arrow::ipc::feather::WriteTable(*table, file.get()));

  ARROW_RETURN_NOT_OK(file->Close());
  return arrow::Status::OK();
}

}  // namespace

bool ConvergenceResultArray::SaveToFile(const std::string& filename) const {
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
