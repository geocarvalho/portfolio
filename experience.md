---
layout: page
title: Experience
permalink: /experience/
---

My professional journey so far.

<div class="experience-list">
{% for job in site.data.experience %}
<div class="experience-item">
  <h3>{{ job.role }}</h3>
  <p class="experience-meta">{{ job.company }} &middot; {{ job.period }}</p>
  <p>{{ job.description }}</p>
</div>
{% endfor %}
</div>
