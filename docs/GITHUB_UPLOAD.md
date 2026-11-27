# 上传到GitHub

仓库内容已整理在 `code/final`。上传本目录的内容即可，不要在包含数GB数据的项目根目录执行 `git add .`。

## 网页上传

1. 在GitHub登录，点击“New repository”，例如命名为 `copd-myeloid`，按投稿要求选择公开。
2. 创建仓库后选择“uploading an existing file”或“Add file → Upload files”。
3. 将本目录中的 `README.md`、`run.py`、requirements文件，以及 `analysis`、`tools`、`resources`、`docs` 上传，保持目录层级。
4. 隐藏文件 `.gitignore` 和 `.github/workflows/check.yml` 也应上传；系统文件管理器需开启“显示隐藏文件”。
5. 提交后检查README中的相对链接，再创建一个Release（例如 `v1.0.0`），供稿件固定版本引用。

如果网页难以处理隐藏文件或多层目录，使用下面的终端方式。

## 终端上传

先在GitHub创建一个**空仓库**，不要自动添加README；把命令中的用户名和仓库名替换为实际值：

```bash
cd /data/sgp/copd/data20261125/code/final
git init -b main
git config user.name "YOUR_GITHUB_USERNAME"
git config user.email "YOUR_VERIFIED_EMAIL_OR_GITHUB_NOREPLY_EMAIL"
git add .
git status --short
git commit -m "Add reproducible COPD myeloid analysis and figures"
git remote add origin https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY.git
git push -u origin main
```

使用GitHub登录授权、SSH密钥或个人访问令牌完成认证；不要把令牌写入代码或remote URL。若本目录已经是Git仓库，省略 `git init`；若origin已存在，先用 `git remote -v` 核对，而不要重复添加。

`user.name` 和 `user.email` 只配置当前仓库。邮箱请使用GitHub已验证邮箱，或GitHub **Settings → Emails** 显示的隐私邮箱。出现 `Author identity unknown` 时，先完成这两项配置，再运行 `git commit` 和 `git push`；`src refspec main does not match any` 通常是前面的提交失败、尚未产生首次提交所致。

官方步骤：[Adding locally hosted code to GitHub](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github)。

## 更新稿件中的代码声明

仓库上传后，用真实公开地址重新生成文稿，避免手动修改后与代码资源不同步：

```bash
python run.py --data-root /data/sgp/copd/data20261125 --stage documents --repository-url https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY
```

该命令重建Word、原始coloc表格及RIS，不重新拟合模型或重绘图。建议同时记录Release标签或提交SHA。投稿前检查：Word中没有“待填写”标记；原始GEO来源完整；所用GWAS、外部coloc汇总及整合对象也有可访问来源。代码的公开URL不能替代缺失的数据来源声明。
