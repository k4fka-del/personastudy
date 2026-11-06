## Kết nối R trên máy với StatNode

StatNode có thể chạy R theo hai cách:

- **Không kết nối (mặc định):** dùng **WebR** — R chạy ngay trong trình duyệt, không cần cài thêm gì.
- **Kết nối với R ở máy:** dùng **R bạn đã cài trên Windows** qua app **StatNode-LocalR** — phù hợp mô hình nặng, package không có trên WebR (ví dụ GARCH, TVAR).

Hướng dẫn dưới đây dành cho chế độ **Kết nối với R ở máy**.


## Bước 1 — Cài R từ CRAN

1. Mở https://cran.r-project.org/bin/windows/base/
2. Tải bản **R mới nhất cho Windows** và cài với tùy chọn mặc định.
3. Kiểm tra trong PowerShell:

```powershell
Rscript --version
```

Nếu lệnh không chạy, mở lại PowerShell sau khi cài hoặc thêm vào PATH thư mục `bin` của R (ví dụ `C:\Program Files\R\R-4.x.x\bin`).


## Bước 2 — Tải và cài StatNode-LocalR

1. [Tải StatNode-LocalR](https://drive.google.com/file/d/1iCJi5AUWMLoGd-xxQ6Q7QmLx0ps_tk4f/view?usp=drive_link).
2. Chạy file installer vừa tải (`StatNode-LocalR_*_x64-setup.exe`).
3. Nếu Windows SmartScreen cảnh báo (chưa ký code): chọn **More info → Run anyway**.
4. Sau khi cài, mở **StatNode-LocalR** từ Start Menu.

Cửa sổ app sẽ hiện **mã ghép nối 6 số**. Giữ app chạy (bấm X chỉ thu nhỏ xuống khay hệ thống, server vẫn hoạt động).


## Bước 3 — Ghép nối trên statnode.io

1. Mở https://statnode.io và đăng nhập (nếu cần).
2. Vào **Cài đặt** (Settings).
3. Mục **Kết nối với R ở máy của bạn** → chọn **Kết nối với R ở máy**.
4. Bấm **Tìm R ở máy**.
5. Nhập **mã 6 số** từ cửa sổ StatNode-LocalR → bấm **Ghép nối**.

Sau khi ghép nối thành công, trang sẽ tải lại và StatNode dùng **R trên máy bạn** thay vì WebR. Bạn không cần chờ đồng bộ package WebR nữa — package dùng thư viện R đã cài trên máy.


## Kiểm tra nhanh

1. Mở tab **R Workspace**, gõ `1 + 1` → kết quả `2`.
2. Chạy workflow có node desktop-only (ví dụ GARCH) — chỉ chạy được khi đã kết nối R ở máy.


## Ngắt kết nối / quay lại WebR

- Trong **Cài đặt** → chọn **Không kết nối** để dùng lại WebR.
- Hoặc bấm **Ngắt kết nối** khi đang ở chế độ R ở máy.
- Trong app StatNode-LocalR có thể bấm **Huỷ ghép nối** để sinh mã mới.


## Lỗi thường gặp

| Triệu chứng | Cách xử lý |
|-------------|------------|
| **Failed to fetch** khi tìm LocalR | Trình duyệt có thể đã chặn quyền mạng cục bộ — xem [Bật lại quyền trình duyệt](/docs/local-r-browser-permission) |
| Không tìm thấy StatNode-LocalR | Mở app từ Start Menu; kiểm tra icon khay hệ thống |
| Không tìm thấy R | Cài R từ CRAN, mở lại StatNode-LocalR |
| Sai mã ghép nối | Mã đổi mỗi 10 phút nếu chưa ghép — copy lại từ app |
| SmartScreen chặn installer | More info → Run anyway (hoặc ký code khi phát hành chính thức) |
