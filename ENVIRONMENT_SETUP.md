# 新电脑环境配置与项目恢复

这份说明用于在新的 Windows 笔记本上恢复本仓库的开发环境。项目采用：

```text
Windows
  -> WSL2 Ubuntu
      -> GCC/G++、CMake、Ninja、GDB、Eigen
  -> CLion
      -> WSL Toolchain
      -> WSL-Debug CMake Profile
```

项目要求：

- CMake 3.20 或更高版本；
- 支持 C++17 的编译器；
- Eigen 3.3 或更高版本（用于 Eigen 示例和 `Project_01`）；
- 推荐使用 Ninja 和 GDB。

构建目录、CLion 本地设置和编译产物都可以重新生成，不需要从旧电脑复制。

---

## 1. 离开旧电脑前

在仓库根目录检查是否还有未提交内容：

```bash
git status
git log -5 --oneline
```

确认最新提交已经推送到 GitHub：

```bash
git push origin main
```

只有本地 commit、没有 push 的内容不能直接在另一台电脑上通过 `git clone` 获得。

不要手工复制下面这些目录：

```text
build/
cmake-build-debug/
cmake-build-wsl-debug/
.idea/
```

它们包含旧电脑的绝对路径、工具链缓存或 IDE 本地设置，在新电脑上重新生成更可靠。

---

## 2. 安装 WSL2 和 Ubuntu

以管理员身份打开 PowerShell，执行：

```powershell
wsl --install -d Ubuntu
```

重启 Windows，首次启动 Ubuntu 时按提示创建 Linux 用户名和密码。输入 Linux 密码时终端不会显示字符，这是正常现象。

在 PowerShell 中确认 Ubuntu 使用 WSL2：

```powershell
wsl --list --verbose
```

预期 `VERSION` 列为 `2`。如果不是，可以执行：

```powershell
wsl --set-version Ubuntu 2
```

官方参考：[Microsoft WSL 安装说明](https://learn.microsoft.com/windows/wsl/install)

---

## 3. 安装 C++ 工具链和 Eigen

打开 Ubuntu/WSL 终端，执行：

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build gdb git libeigen3-dev
```

检查工具是否可用：

```bash
g++ --version
cmake --version
ninja --version
gdb --version
git --version
```

检查 Eigen 的 CMake 配置文件：

```bash
test -f /usr/share/eigen3/cmake/Eigen3Config.cmake && echo "Eigen is ready"
```

如果输出 `Eigen is ready`，说明依赖已安装。

---

## 4. 配置 Git 和 GitHub SSH

在 WSL 中设置提交身份。把示例内容替换为自己的姓名和 GitHub 邮箱：

```bash
git config --global user.name "Your Name"
git config --global user.email "your-email@example.com"
```

查看是否已有 SSH 密钥：

```bash
ls -la ~/.ssh
```

新电脑通常需要生成新密钥：

```bash
ssh-keygen -t ed25519 -C "your-email@example.com"
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
cat ~/.ssh/id_ed25519.pub
```

只把 `id_ed25519.pub` 的内容添加到 GitHub 的 SSH keys 页面。不要复制或分享私钥 `id_ed25519`。

测试连接：

```bash
ssh -T git@github.com
```

官方参考：[GitHub 生成并添加 SSH 密钥](https://docs.github.com/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent?platform=linux)

---

## 5. 克隆仓库

为了延续当前工作方式，可以把仓库放在 Windows 的 `C:\code`，它在 WSL 中对应 `/mnt/c/code`：

```bash
mkdir -p /mnt/c/code
cd /mnt/c/code
git clone git@github.com:47-5/learning_cpp_for_chemistry.git
cd learning_cpp_for_chemistry
```

检查分支和状态：

```bash
git branch --show-current
git status
```

预期分支为 `main`，工作区没有未提交改动。

如果更重视 WSL 文件系统性能，也可以克隆到 `~/code`。无论选择哪个位置，都不要在后续移动已经生成过的 CMake 构建目录；移动源码后重新配置即可。

---

## 6. 先用命令行验证项目

在 WSL 的仓库根目录执行：

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

其中：

- `-S .` 指定源码目录；
- `-B build` 指定可重新生成的构建目录；
- `-G Ninja` 使用 Ninja；
- `CMAKE_BUILD_TYPE=Debug` 生成适合调试的程序；
- `ctest` 运行已注册的测试。

构建成功后，可以运行：

```bash
./build/00_cpp_fundamentals/study_hello_world
./build/00_cpp_fundamentals/fizzbuzz
./build/01_basic/basic_hello
./build/02_eigen/eigen_matrix_demo
./build/03_miniHF/minihf_placeholder
```

`project01_molecule` 是静态库目标，因此没有可直接运行的同名程序，但它应当在完整构建时成功编译。

---

## 7. 配置 CLion 的 WSL Toolchain

安装并启动 CLion，然后打开 Windows 路径：

```text
C:\code\learning_cpp_for_chemistry
```

进入：

```text
Settings
  -> Build, Execution, Deployment
  -> Toolchains
```

添加 `WSL` Toolchain，选择 Ubuntu。确认 CLion 能检测到：

```text
CMake        /usr/bin/cmake
Build Tool   /usr/bin/ninja
C Compiler   /usr/bin/cc
C++ Compiler /usr/bin/c++
Debugger     /usr/bin/gdb
```

官方参考：[CLion WSL2 Toolchain](https://www.jetbrains.com/help/clion/how-to-use-wsl-development-environment-in-product.html)

---

## 8. 配置 CLion 的 CMake Profile

进入：

```text
Settings
  -> Build, Execution, Deployment
  -> CMake
```

创建或保留下面的 Profile：

```text
Name:             WSL-Debug
Enable profile:   Yes
Build type:       Debug
Toolchain:        WSL
Generator:        Ninja
Build directory:  cmake-build-wsl-debug
```

禁用不使用的 Windows `Debug` Profile，否则 CLion 会同时生成：

```text
cmake-build-debug/
cmake-build-wsl-debug/
```

两个目录并不冲突，但只使用 WSL 时没有必要同时保留。CLion 官方也建议禁用当前不用的 Profile，以减少重新加载时间和潜在错误。

配置完成后执行：

```text
Tools -> CMake -> Reset Cache and Reload Project
```

等待输出出现：

```text
Build files have been written to: .../cmake-build-wsl-debug
```

然后选择 `study_hello_world` 或其他 target 进行运行和调试。

官方参考：[CLion CMake Profiles](https://www.jetbrains.com/help/clion/cmake-profile.html)

---

## 9. 日常工作流程

开始工作前：

```bash
git pull --ff-only
git status
```

完成一个小阶段后：

```bash
git status
git diff
git add path/to/file
git commit -m "描述本次学习内容"
git push origin main
```

不要提交 API token、密码、SSH 私钥、CLion 本地缓存或 CMake 构建目录。

---

## 10. 常见问题

### 删除构建目录后 CLion 报 `No such file or directory`

构建目录可以删除，但删除后需要先重新生成：

```text
Tools -> CMake -> Reset Cache and Reload Project
```

不要直接点击 Build，因为此时还没有 Ninja 构建文件。

### 重新加载后出现两个构建目录

说明 Windows `Debug` 和 `WSL-Debug` 两个 CMake Profile 都处于启用状态。在 CMake Settings 中禁用不需要的 Profile。

### CMake 提示没有找到 Eigen

在 WSL 中执行：

```bash
sudo apt update
sudo apt install libeigen3-dev
```

然后 Reset Cache and Reload Project。

### 从旧电脑复制项目后出现旧路径

删除明确的生成目录，例如：

```bash
rm -rf build
```

然后重新执行 CMake 配置。不要删除仓库根目录，也不要复制旧电脑的 `CMakeCache.txt`。

### GitHub SSH 返回 `Permission denied (publickey)`

检查：

```bash
ssh-add -l
ssh -T git@github.com
```

确认新电脑的公钥已经添加到 GitHub，并且使用的是 WSL 内的 `~/.ssh`。

---

## 11. 换机完成检查表

- [ ] WSL2 Ubuntu 可以启动；
- [ ] `g++`、`cmake`、`ninja`、`gdb`、`git` 可用；
- [ ] Eigen 配置文件存在；
- [ ] GitHub SSH 测试成功；
- [ ] 仓库已从 GitHub 克隆；
- [ ] 命令行完整构建成功；
- [ ] CTest 全部通过；
- [ ] CLion WSL Toolchain 检测正常；
- [ ] 仅启用需要的 CMake Profile；
- [ ] CLion 可以运行并调试 `study_hello_world`。
