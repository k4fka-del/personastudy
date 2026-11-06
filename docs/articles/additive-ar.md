# Additive AR (aar)

## 1. Mô hình / công cụ này dùng để làm gì?

Mô hình additive AR ước lượng phi tuyến bằng tổng các hàm một chiều của từng lag, tránh giả định dạng tham số cụ thể như SETAR/LSTAR.

Node này yêu cầu StatNode Desktop vì package tsDyn chưa có bản WebR.

## 2. Nền tảng lý thuyết


$$y_t = \mu + \sum_{j=1}^{p} f_j(y_{t-j}) + \varepsilon_t$$

các $f_j$ được ước lượng phi tham số (splines) trong tsDyn::aar().

## 3. Điều kiện sử dụng


- Phù hợp khám phá phi tuyến khi chưa có giả thuyết regime rõ.
- So sánh AIC/BIC với AR tuyến tính để đánh giá cải thiện.

## 4. Quy trình thực hiện (Workflow)


1. aar() với embedding $m$.
2. Chẩn đoán phần dư.
3. So sánh với linear().
4. Dự báo.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình Additive AR (aar)](/docs/images/additive-ar-settings.png)

*Hình 1. Panel Settings của node `additive-ar` trên canvas StatNode (Desktop).*

### Output mẫu

![Output Additive AR (aar)](/docs/images/additive-ar-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
=== aar-main ===
--- Additive AR Summary ---
Embedding dimension m = 3

Smooth functions f_j estimated for lags 1–3.
Residual SE: 0.412
AIC: 415.2   BIC: 431.8

Ljung-Box on residuals (lag 12): p = 0.34 (no autocorrelation)

Variable: lynx
```


## 6. Diễn giải kết quả

Nếu AAR cải thiện tiêu chí thông tin so với AR tuyến tính, có bằng chứng phi tuyến trong động lực lag.

## 7. Tài liệu tham khảo


- Chen, R., & Tsay, R. S. (1993). Nonlinear additive AR models. *JASA*.
- Stigler, M., et al. tsDyn vignette "Nonlinear time series in practice".
