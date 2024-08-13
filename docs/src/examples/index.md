```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const tutorials = [
  {
    href: "generated/introduction",
    src: "../introduction.png",
    caption: "Introduction",
    desc: "Groups, elements, representations",
  },
];
</script>
```

# Examples

## Tutorials

```@raw html
<Gallery :images="tutorials" />
```
