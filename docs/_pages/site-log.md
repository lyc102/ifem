---
permalink: /logs/site/
title: "iFEM site development notes"
sidebar:
    nav: docs
---



## Adding new pages in documentation site
- Writing in markdown.
- To put it under a certain category, copy the markdown file together with its image folder (under the same name) to `_pages` (documents for iFEM usage, community research, update, logs, etc). Directly copying `_todo` folders `fem`, `afem`, `solver` to `_fem`, `_afem`, `_solver` is easier but we need to add those folders to `include:` in `_config.yml`. 
- Adding the following preamble to the markdown file: where `category` is either `fem`, `afem`, `solver`, etc.
    ```html
    ---
    permalink: /category/page-name/
    title: "Page Title"
    sidebar:
        nav: docs
    ---
    ```
- The default left sidebar navigation is `docs`, there are other navigation links in `_data/navigation.yml`.

## Testing pages
- To test building the site locally, install `ruby` and `jekyll` (see [below](#july-25-2021) for details).
- To add a test page, we can add `test_something.md` in `_pages`, then the page can be accessed directly without being shown up in the site using the `permalink` in the markdown file's preamble. If `permalink: /test/test_something/` is added, then `http://localhost:4000/ifem/test/test_something/` should work in the local served version and `http://lyc1012.github.io/ifem/test/test_something/` when deployed.


## CHANGELOG

### July 27, 2021
- [x] Adding Jekyll preambles to un-deployed `.md` converted from the old `/ifemdoc/*/*.ipynb` are in `_todo`.
- [x] re-organizing `_todo` folders to their corresponding iFEM main repo folders.
- [x] manually checked the links in nav.

### July 25, 2021
Manually verified building on MacOS. Installation guide on MacOS: [https://jekyllrb.com/docs/installation/macos/](https://jekyllrb.com/docs/installation/macos/). The easiest way is to use `brew` to install `ruby` and `jekyll`.
```bash
brew install ruby
```
then
```bash
# If you're using Zsh
echo 'export PATH="/usr/local/opt/ruby/bin:/usr/local/lib/ruby/gems/3.0.0/bin:$PATH"' >> ~/.zshrc

# If you're using Bash
echo 'export PATH="/usr/local/opt/ruby/bin:/usr/local/lib/ruby/gems/3.0.0/bin:$PATH"' >> ~/.bash_profile
```
After Ruby is installed, install Jekyll using `gem install --user-install bundler jekyll`. Now under the `docs` folder, run `bundle install` to install the dependencies then `bundle exec jekyll serve` will start a localhost to serve the pages at `http://127.0.0.1:4000/ifem/`.

### July 23, 2021
- Set up the doc site.
- Navigation `/_data/navigation.yml` controls all the navigation bars.
- Main scss is `assets/css/main.scss` controls the main CSS styles. Global var is defined before `@import`, syntax highlighting and colors can be defined in the same file after `@import "minimal-mistakes/skins/{{ site.minimal_mistakes_skin | default: 'default' }}"`, adjustment to others or customized CSS can be added after `@import "minimal-mistakes"`.
- $\LaTeX$ can be added using the usual `$ $` or `$$ $$`. 
    - Kramdown supports `\eqref` as well. For an example, please refer to `_docs/fem-poisson.md`. Inside the mathmode, at the end of an equation, add `\tag{num}\label{label_name}`, then `\eqref{}` can be directly called in the main text.
    - In MathJax, `\{\}` needs to be rewritten to `\\{\\}` to escape the backslash interpreter for HTML; also for nested subscripts or subscripts in superscripts, the underscore `_` needs a backslash `\` before it, for example, `$u_{g_D}$` needs to be `$u_{g\_D}$`.