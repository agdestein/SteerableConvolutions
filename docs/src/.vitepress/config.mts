import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import { transformerMetaWordHighlight } from '@shikijs/transformers';

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  // head: [['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }]],

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
        md.use(mathjax3),
        md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
    codeTransformers: [transformerMetaWordHighlight(),],
  },

  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    // logo: { src: '/logo.svg' },
    logo: {
      'light': '/logo.svg',
      'dark': '/logo.svg'
    },
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and <a href="https://luxdl.github.io/DocumenterVitepress.jl/stable/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    },
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Getting Started', link: '/getting_started' },
      {
        text: 'Examples',
        items: [
          { text: 'Examples', link: '/examples/' },
          {
            text: 'Tutorials',
            items: [
              { text: 'Groups and representations', link: '/examples/generated/introduction' },
              { text: 'Equivariance', link: '/examples/generated/equivariance' },
              // { text: 'Scratch', link: '/examples/generated/scratch' },
            ],
          },
        ],
      },
      { text: 'API reference', link: '/api' },
      {
        text: 'About',
        items: [
          { text: 'About', link: '/about/' },
          { text: 'License', link: '/about/license' },
          { text: 'Citing', link: '/about/citing' },
          { text: 'Contributing', link: '/about/contributing' },
          { text: 'Package versions', link: '/about/versions' },
        ],
      },
      { text: 'References', link: '/references' },
    ],
    sidebar: {
      "/examples/": {
        text: 'Examples',
        items: [
          { text: 'Examples', link: '/examples/' },
          {
            text: 'Tutorials',
            items: [
              { text: 'Groups and representations', link: '/examples/generated/introduction' },
              { text: 'Equivariance', link: '/examples/generated/equivariance' },
              // { text: 'Scratch', link: '/examples/generated/scratch' },
            ],
          },
        ],
      },
      "/about/": {
        text: 'About',
        items: [
          { text: 'About', link: '/about/' },
          { text: 'License', link: '/about/license' },
          { text: 'Citing', link: '/about/citing' },
          { text: 'Contributing', link: '/about/contributing' },
          { text: 'Package versions', link: '/about/versions' },
        ],
      },
    },
  },
})
