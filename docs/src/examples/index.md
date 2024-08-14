```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const tutorials = [
  {
    href: "generated/introduction",
    src: "../introduction.png",
    caption: "Groups and representations",
    desc: "Introductory tutorial to groups and representations",
  },
  {
    href: "generated/equivariance",
    src: "../introduction.png",
    caption: "Equivariance",
    desc: "Introduction to geometric transforms and equivariance",
  },
];
</script>
```

# Examples

## Tutorials

```@raw html
<Gallery :images="tutorials" />
```
