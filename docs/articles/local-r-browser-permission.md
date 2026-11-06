## Bật lại quyền trình duyệt cho StatNode-LocalR

Khi statnode.io kết nối với **StatNode-LocalR** trên máy bạn, trình duyệt cần quyền truy cập **mạng cục bộ / loopback** (localhost). Nếu lỡ bấm **Chặn** ở hộp thoại *“Truy cập vào các ứng dụng và dịch vụ khác trên thiết bị này”*, StatNode sẽ báo lỗi **Failed to fetch** và không tìm thấy LocalR dù app đang chạy.


## Chrome hoặc Edge (Windows)

1. Mở https://statnode.io và vào **Cài đặt** → **Kết nối với R ở máy**.
2. Bấm icon **cài đặt site** (hai thanh ngang có chấm tròn) **bên trái** thanh địa chỉ, cạnh `statnode.io`.
3. Tìm mục **Truy cập mạng cục bộ** / **Local network access** (hoặc tương tự).
4. Đổi từ **Chặn** → **Cho phép**.
5. **Tải lại** trang (F5 hoặc Ctrl+R).
6. Bấm lại **Tìm R ở máy** và ghép nối mã 6 số từ StatNode-LocalR.

Nếu không thấy mục trên, mở **Cài đặt trang web** từ menu đó → tìm **Quyền** → **Truy cập mạng cục bộ** → **Cho phép** cho statnode.io.


## Firefox

1. Bấm icon **ổ khóa / thông tin trang** bên trái URL.
2. Mở **Thông tin trang** hoặc **Quyền**.
3. Tìm quyền liên quan **mạng cục bộ** / **local network** và đặt **Cho phép**.
4. Tải lại trang và thử kết nối lại.


## Kiểm tra nhanh

- StatNode-LocalR đang mở (icon khay hệ thống hoặc cửa sổ app).
- Trong **Cài đặt**, chọn **Kết nối với R ở máy** → **Tìm R ở máy** không còn báo **Failed to fetch**.
- Ghép nối thành công → card **Kết nối với R ở máy** hiện *Đang sử dụng R ở máy*.


## Vẫn không được?

| Triệu chứng | Cách xử lý |
|-------------|------------|
| Vẫn Failed to fetch sau khi Cho phép | Tải lại cứng (Ctrl+Shift+R); đóng tab statnode.io mở lại |
| Không thấy mục quyền mạng cục bộ | Cập nhật Chrome/Edge lên bản mới; thử trình duyệt khác |
| Không tìm thấy LocalR | Xem [Hướng dẫn kết nối R trên máy](/docs/local-r) — cài và mở StatNode-LocalR trước |
