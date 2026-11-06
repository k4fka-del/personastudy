# LSTAR (Logistic Smooth Transition Autoregressive Model)

## 1. Mô hình / công cụ này dùng để làm gì?

LSTAR mô tả chuyển đổi mượt giữa hai hành vi tự hồi quy thay vì nhảy cóc như SETAR.

Node này yêu cầu StatNode Desktop vì package tsDyn chưa có bản WebR.

## 2. Nền tảng lý thuyết


$$y_t = (\phi_0' x_t)(1-G(z_{t-d};\gamma,c)) + (\psi_0' x_t) G(z_{t-d};\gamma,c) + \varepsilon_t$$

$$G(z;\gamma,c) = \left(1+\exp(-\gamma(z-c))\right)^{-1}$$

## 3. Điều kiện sử dụng


- Cần mẫu đủ dài; tham số $\gamma$ lớn có thể gần SETAR.
- Kiểm định LM (star) để xác nhận phi tuyến chuyển tiếp mượt.

## 4. Quy trình thực hiện (Workflow)


1. Ước lượng AR tuyến tính tham chiếu.
2. star() kiểm định số regime / tuyến tính.
3. lstar() ước lượng.
4. Trực quan hóa hàm chuyển tiếp và dự báo.

## 5. Thực hành trên phần mềm StatNode

### Cấu hình node

![Cấu hình LSTAR (Logistic Smooth Transition Autoregressive Model)](/docs/images/lstar-settings.png)

*Hình 1. Panel Settings của node `lstar` trên canvas StatNode (Desktop).*

### Output mẫu

![Output LSTAR (Logistic Smooth Transition Autoregressive Model)](/docs/images/lstar-output.png)

*Hình 2. Output sau khi chạy workflow.*


```
=== lstar-main ===
--- LSTAR Summary ---
Logistic Smooth Transition AR
Gamma: 12.4 | Threshold c: 7.38

Low regime AR coefficients:
  L1: 0.598***  L2: -0.091

High regime AR coefficients:
  L1: 0.521***  L2: 0.142*

Residual variance: 0.178
AIC: 408.9

Variable: lynx
```


## 6. Diễn giải kết quả

Tham số $\gamma$ đo độ dốc chuyển tiếp; $c$ là điểm giữa logistic. Regime mờ hơn SETAR.

## 7. Tài liệu tham khảo


- Teräsvirta, T. (1994). Specification, estimation, and evaluation of smooth transition autoregressive models. *Journal of the American Statistical Association*.
- van Dijk, D., Teräsvirta, T., & Franses, P. H. (2002). Smooth transition autoregressive models. *Econometric Reviews*.
