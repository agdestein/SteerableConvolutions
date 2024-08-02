```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const examples = [
  {
    href: "generated/example",
    src: "../example.png",
    caption: "Example",
    desc: "Steerable conv example",
  },
];
</script>
```

# Examples

```@raw html
<Gallery :images="examples" />
```
