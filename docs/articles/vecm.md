# VECM (Vector Error Correction Model)

## 1. Mô hình / công cụ này dùng để làm gì?

VECM mô tả động lực ngắn hạn của hệ nhiều biến chuỗi thời gian có quan hệ cointegration, đồng thời áp đặt ràng buộc điều chỉnh sai số dài hạn.

Node này yêu cầu StatNode Desktop vì package tsDyn chưa có bản WebR.

## 2. Nền tảng lý thuyết


Hệ VECM($p$) với rank $r$:

$$\Delta Y_t = \Pi Y_{t-1} + \sum_{i=1}^{p-1} \Gamma_i \Delta Y_{t-i} + \Phi D_t + \varepsilon_t$$

trong đó $\Pi = \alpha \beta'$, $\beta$ chứa các vector cointegration, $\alpha$ là hệ số điều chỉnh tốc độ.

## 3. Điều kiện sử dụng


- Các biến nên là I(1) và cointegrated (kiểm tra Johansen trước).
- Chọn độ trễ $p$ và rank $r$ dựa trên thông tin mẫu và kiểm định.
- Mẫu đủ dài so với số tham số ước lượng.

## 4. Quy trình thực hiện (Workflow)


1. Kiểm tra tính dừng / tích phân bậc một.
2. Johansen xác định rank cointegration.
3. Ước lượng VECM (ML hoặc 2OLS) với tsDyn::VECM().
4. Kiểm tra ổn định, nhân quả Granger, IRF và FEVD.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình VECM (Vector Error Correction Model)](/docs/images/vecm-settings.png)

*Hình 1. Panel Settings của node `vecm` trên canvas StatNode (Desktop).*

### Output mẫu

![Output VECM (Vector Error Correction Model)](/docs/images/vecm-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
=== vecm-main ===
--- VECM Summary ---
Vector Error Correction Model (VECM)
Lag order: 2
Cointegration rank: 1
Estimator: ML (Johansen)

Adjustment coefficients (alpha):
              long.run short.run
long.run     -0.042***    0.011
short.run     0.018      -0.035**

Cointegrating vector (beta):
long.run = 1.000
short.run = -0.982***

Short-run dynamics and deterministic terms reported in full summary.

Variables: long.run, short.run
Cointegration rank used: 1
```


## 6. Diễn giải kết quả

Hệ số điều chỉnh âm và có ý nghĩa cho thấy biến hồi phục về quan hệ dài hạn; hệ số ngắn hạn mô tả phản ứng tức thời.

## 7. Tài liệu tham khảo


- Engle, R. F., & Granger, C. W. J. (1987). Co-integration and error correction. *Econometrica*.
- Johansen, S. (1991). Estimation and hypothesis testing of cointegration vectors. *Econometrica*.
- Lütkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer.
- Stigler, M. (2024). tsDyn: Nonlinear Time Series Models with Regime Switching. CRAN.
