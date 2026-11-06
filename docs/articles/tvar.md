# TVAR (Threshold Vector Autoregression)

## 1. Mô hình / công cụ này dùng để làm gì?

TVAR mở rộng VAR bằng cách cho phép hệ số thay đổi theo regime đa biến, hữu ích khi chính sách hoặc shock lớn tạo hai trạng thái hệ thống.

Node này yêu cầu StatNode Desktop vì package tsDyn chưa có bản WebR.

## 2. Nền tảng lý thuyết


Với $K$ biến và 2 regime:

$$Y_t = A_{1,0} + \sum_{i=1}^{p} A_{1,i} Y_{t-i} + \varepsilon_t \quad \text{nếu } s_{t-d} \leq \gamma$$

$$Y_t = A_{2,0} + \sum_{i=1}^{p} A_{2,i} Y_{t-i} + \varepsilon_t \quad \text{nếu } s_{t-d} > \gamma$$

## 3. Điều kiện sử dụng


- Ít nhất hai biến; mẫu dài.
- Grid search chậm — giới hạn nthresh và trim hợp lý.
- TVAR.LRtest kiểm định tuyến tính trước.

## 4. Quy trình thực hiện (Workflow)


1. VAR tuyến tính tham chiếu.
2. TVAR.LRtest.
3. TVAR() ước lượng.
4. IRF theo regime, phân loại regime, dự báo.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình TVAR (Threshold Vector Autoregression)](/docs/images/tvar-settings.png)

*Hình 1. Panel Settings của node `tvar` trên canvas StatNode (Desktop).*

### Output mẫu

![Output TVAR (Threshold Vector Autoregression)](/docs/images/tvar-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
=== tvar-main ===
--- TVAR Summary ---
Threshold Vector Autoregressive Model
Lag: 2 | Regimes: 2
Threshold gamma: -0.842

Regime 1 coefficients (excerpt):
  long.run L1: 0.451***
  short.run L1: 0.112

Regime 2 coefficients (excerpt):
  long.run L1: 0.389**
  short.run L1: 0.204*

Log-likelihood: -1842.6
AIC: 3721.2

Variables: long.run, short.run
```


## 6. Diễn giải kết quả

Hệ số khác biệt giữa regime cho thấy truyền dẫn shock phụ thuộc trạng thái hệ thống.

## 7. Tài liệu tham khảo


- Tsay, R. S. (1998). Testing and modeling multivariate threshold models. *Journal of the American Statistical Association*.
- Balke, N. S. (2000). Credit and economic activity. *Federal Reserve Bank of Dallas Economic Review*.
