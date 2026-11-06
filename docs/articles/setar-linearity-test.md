# SETAR • Linearity Test

## 1. Mô hình / công cụ này dùng để làm gì?

setarTest kiểm định AR tuyến tính vs SETAR.

Sub-node này gắn với node cha **setar** trên StatNode Desktop (package tsDyn).

## 2. Nền tảng lý thuyết

Hansen (1996) bootstrap LM-type test.

## 3. Điều kiện sử dụng

- Dữ liệu chuỗi thời gian đủ dài, ít missing.
- Node chỉ chạy trên StatNode Desktop.
- Nên hoàn thành các bước tiền đề (kiểm định, chọn lag) trước sub-node.

## 4. Quy trình thực hiện (Workflow)


1. Kéo node cha (setar) và Input Data trên canvas StatNode Desktop.
2. Thêm sub-node **SETAR • Linearity Test** từ menu sub-node.
3. Cấu hình tham số (nếu có) và chạy workflow.
4. Đọc output trong panel kết quả.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình SETAR • Linearity Test](/docs/images/setar-linearity-test-settings.png)

*Hình 1. Panel Settings của node `setar-linearity-test` trên canvas StatNode (Desktop).*

### Output mẫu

![Output SETAR • Linearity Test](/docs/images/setar-linearity-test-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
=== setar-linearity-test ===
--- SETAR Linearity Test (setarTest) ---
Test: 1 vs SETAR (1 threshold)
Bootstrap replications: 50

Test statistic: 18.42
p-value: 0.024

Conclusion: Reject linear AR at 5% — evidence for threshold nonlinearity.
```


## 6. Diễn giải kết quả

Diễn giải output theo ngữ cảnh mô hình cha: kiểm tra ý nghĩa thống kê, dấu hệ số và so sánh với lý thuyết trước khi kết luận chính sách.

## 7. Tài liệu tham khảo


- Stigler, M. (2024). *tsDyn* package documentation. CRAN.
- Lütkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer.
