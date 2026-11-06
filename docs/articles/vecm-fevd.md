# VECM • FEVD

## 1. Mô hình / công cụ này dùng để làm gì?

Phân rã phương sai sai số dự báo của một biến do shock từng biến trong hệ.

Sub-node này gắn với node cha **vecm** trên StatNode Desktop (package tsDyn).

## 2. Nền tảng lý thuyết

FEVD$_j(h)$ = phần trăm MSE horizon $h$ của biến $i$ do shock biến $j$.

## 3. Điều kiện sử dụng

- Dữ liệu chuỗi thời gian đủ dài, ít missing.
- Node chỉ chạy trên StatNode Desktop.
- Nên hoàn thành các bước tiền đề (kiểm định, chọn lag) trước sub-node.

## 4. Quy trình thực hiện (Workflow)


1. Kéo node cha (vecm) và Input Data trên canvas StatNode Desktop.
2. Thêm sub-node **VECM • FEVD** từ menu sub-node.
3. Cấu hình tham số (nếu có) và chạy workflow.
4. Đọc output trong panel kết quả.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình VECM • FEVD](/docs/images/vecm-fevd-settings.png)

*Hình 1. Panel Settings của node `vecm-fevd` trên canvas StatNode (Desktop).*

### Output mẫu

![Output VECM • FEVD](/docs/images/vecm-fevd-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
(Kết nối sub-node sau node cha trên canvas Desktop và chạy workflow để xem output.)
```


## 6. Diễn giải kết quả

Diễn giải output theo ngữ cảnh mô hình cha: kiểm tra ý nghĩa thống kê, dấu hệ số và so sánh với lý thuyết trước khi kết luận chính sách.

## 7. Tài liệu tham khảo


- Stigler, M. (2024). *tsDyn* package documentation. CRAN.
- Lütkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer.
