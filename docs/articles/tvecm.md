# TVECM (Threshold Vector Error Correction Model)

## 1. Mô hình / công cụ này dùng để làm gì?

TVECM kết hợp cointegration với ngưỡng: tốc độ điều chỉnh về cân bằng dài hạn có thể khác nhau giữa regime.

Node này yêu cầu StatNode Desktop vì package tsDyn chưa có bản WebR.

## 2. Nền tảng lý thuyết


Mở rộng VECM với $\Pi$ phụ thuộc regime $s_t$:

$$\Delta Y_t = \Pi^{(j)} Y_{t-1} + \sum_{i=1}^{p-1} \Gamma_i^{(j)} \Delta Y_{t-i} + \varepsilon_t \quad \text{khi } s_t = j$$

## 3. Điều kiện sử dụng


- Biến I(1) và có cointegration (kiểm tra Johansen / Seo).
- TVECM.SeoTest và TVECM.HStest cho giả thuyết ngưỡng trên cointegration.

## 4. Quy trình thực hiện (Workflow)


1. VECM tuyến tính.
2. Seo / Hansen–Seo tests.
3. TVECM() ước lượng.
4. IRF và regime.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình TVECM (Threshold Vector Error Correction Model)](/docs/images/tvecm-settings.png)

*Hình 1. Panel Settings của node `tvecm` trên canvas StatNode (Desktop).*

### Output mẫu

![Output TVECM (Threshold Vector Error Correction Model)](/docs/images/tvecm-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
=== tvecm-main ===
--- TVECM Summary ---
Threshold Vector Error Correction Model
Lag: 1 | Rank r: 1 | Regimes: 2
Threshold: -0.91

Regime 1 adjustment (alpha):
  long.run: -0.051***
  short.run: 0.009

Regime 2 adjustment (alpha):
  long.run: -0.028*
  short.run: -0.041**

Cointegrating relationship stable across regimes.
Log-likelihood: -1756.4

Variables: long.run, short.run
```


## 6. Diễn giải kết quả

Reject linear cointegration ủng hộ mô hình có ngưỡng; hệ số $\alpha$ theo regime cho biết tốc độ điều chỉnh không đồng nhất.

## 7. Tài liệu tham khảo


- Hansen, B. E., & Seo, B. (2002). Testing for two-regime threshold cointegration. *Journal of Econometrics*.
- Seo, M. (2006). Bootstrap testing for the null of no cointegration in a threshold vector error correction model. *Journal of Econometrics*.
