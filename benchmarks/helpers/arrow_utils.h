#include <arrow/api.h>

#include <vector>

namespace bench {

inline arrow::Result<std::shared_ptr<arrow::Array>> MakeDoubleArray(
    const std::vector<double>& data) {
  arrow::DoubleBuilder builder;
  ARROW_RETURN_NOT_OK(builder.AppendValues(data.data(), data.size()));
  std::shared_ptr<arrow::Array> array;
  ARROW_RETURN_NOT_OK(builder.Finish(&array));
  return array;
}

// Assumes int maps to arrow::int32().
inline arrow::Result<std::shared_ptr<arrow::Array>> MakeInt32Array(
    const std::vector<int>& data) {
  arrow::Int32Builder builder;
  ARROW_RETURN_NOT_OK(builder.AppendValues(data.data(), data.size()));
  std::shared_ptr<arrow::Array> array;
  ARROW_RETURN_NOT_OK(builder.Finish(&array));
  return array;
}

// Specialization for std::vector<bool> because of its packed nature.
inline arrow::Result<std::shared_ptr<arrow::Array>> MakeBooleanArray(
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

}  // namespace bench
