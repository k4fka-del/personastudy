# TVECM • Hansen–Seo Test

## 1. Mô hình / công cụ này dùng để làm gì?

TVECM.HStest sup-LR với bootstrap.

Sub-node này gắn với node cha **tvecm** trên StatNode Desktop (package tsDyn).

## 2. Nền tảng lý thuyết

Xem tài liệu tsDyn và các paper được trích dẫn ở mục 7.

## 3. Điều kiện sử dụng

- Dữ liệu chuỗi thời gian đủ dài, ít missing.
- Node chỉ chạy trên StatNode Desktop.
- Nên hoàn thành các bước tiền đề (kiểm định, chọn lag) trước sub-node.

## 4. Quy trình thực hiện (Workflow)


1. Kéo node cha (tvecm) và Input Data trên canvas StatNode Desktop.
2. Thêm sub-node **TVECM • Hansen–Seo Test** từ menu sub-node.
3. Cấu hình tham số (nếu có) và chạy workflow.
4. Đọc output trong panel kết quả.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình TVECM • Hansen–Seo Test](/docs/images/tvecm-hs-test-settings.png)

*Hình 1. Panel Settings của node `tvecm-hs-test` trên canvas StatNode (Desktop).*

### Output mẫu

![Output TVECM • Hansen–Seo Test](/docs/images/tvecm-hs-test-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
(Kết nối sub-node sau node cha trên canvas Desktop và chạy workflow để xem output.)
```


## 6. Diễn giải kết quả

Diễn giải output theo ngữ cảnh mô hình cha: kiểm tra ý nghĩa thống kê, dấu hệ số và so sánh với lý thuyết trước khi kết luận chính sách.

## 7. Tài liệu tham khảo


- Stigler, M. (2024). *tsDyn* package documentation. CRAN.
- Lütkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer.
