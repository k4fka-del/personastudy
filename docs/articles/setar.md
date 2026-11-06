# SETAR (Self-Exciting Threshold Autoregressive Model)

## 1. Mô hình / công cụ này dùng để làm gì?

SETAR cho phép hệ số tự hồi quy thay đổi theo regime khi biến ngưỡng vượt qua ngưỡng, phù hợp chuỗi có bất đối xứng và hiệu ứng phi tuyến.

Node này yêu cầu StatNode Desktop vì package tsDyn chưa có bản WebR.

## 2. Nền tảng lý thuyết


Mô hình SETAR(2) regime:

$$y_t = \begin{cases}
\phi_{1,0} + \sum_{j=1}^{p_1} \phi_{1,j} y_{t-j} + \varepsilon_t & z_{t-d} \leq \gamma \\
\phi_{2,0} + \sum_{j=1}^{p_2} \phi_{2,j} y_{t-j} + \varepsilon_t & z_{t-d} > \gamma
\end{cases}$$

## 3. Điều kiện sử dụng


- Chuỗi đủ dài để ước lượng hai (hoặc ba) regime.
- Nên kiểm định tuyến tính (setarTest) trước khi diễn giải SETAR.
- Lựa chọn độ trễ ngưỡng và biến ngưỡng cần cơ sở kinh tế/lý thuyết.

## 4. Quy trình thực hiện (Workflow)


1. Kiểm định tuyến tính (sub-node setar-linearity-test).
2. Chọn đặc tả lag/order (selectSETAR nếu cần).
3. Ước lượng setar().
4. Phân loại regime, dự báo bootstrap.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình SETAR (Self-Exciting Threshold Autoregressive Model)](/docs/images/setar-settings.png)

*Hình 1. Panel Settings của node `setar` trên canvas StatNode (Desktop).*

### Output mẫu

![Output SETAR (Self-Exciting Threshold Autoregressive Model)](/docs/images/setar-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
=== setar-main ===
--- SETAR Summary ---
Self Threshold Autoregressive Model (SETAR)
Number of regimes: 2
Threshold: 7.412

Regime 1 (Low):
  Intercept: 2.145
  AR(1): 0.612***
  AR(2): -0.084
  sigma: 0.421

Regime 2 (High):
  Intercept: 1.892
  AR(1): 0.548***
  AR(2): 0.127*
  sigma: 0.389

AIC: 412.3   BIC: 438.7

Variable: lynx
```


## 6. Diễn giải kết quả

Ngưỡng $\gamma$ chia mẫu thành regime; hệ số khác biệt giữa regime cho thấy phản ứng không đối xứng.

## 7. Tài liệu tham khảo


- Tong, H. (1990). *Non-linear Time Series: A Dynamical System Approach*. Oxford.
- Hansen, B. E. (1999). Testing for linearness. *Journal of Economic Surveys*.
- Chen, R., & Tsay, R. S. (1993). Nonlinear additive AR models. *Journal of the American Statistical Association*.
