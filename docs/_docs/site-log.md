---
permalink: /docs/site-logs/
title: "iFEM site development notes"
sidebar:
    nav: docs
---

## July 23, 2021
- Set up the doc site.
- [Navigation](../_data/navigation.yml) controls the navigation bar.
- [Main scss](../assets/css/main.scss) controls the main CSS styles. Global var is defined before `@import`, syntax highlighting and colors can be defined in the same file after `@import "minimal-mistakes/skins/{{ site.minimal_mistakes_skin | default: 'default' }}"`, adjustment to others or customized CSS can be added after `@import "minimal-mistakes"`.
- In MathJax, `\{\}` needs to be rewritten to `\\{\\}` to escape the backslash interpreter for HTML; also for nested subscripts or subscripts in superscripts, the underscore `_` needs a backslash `\` before it, for example, `$u_{g_D}$` needs to be `$u_{g\_D}$`.