---
layout: page
title: Blog
permalink: /blog/
---

Thoughts, tutorials, and notes on software development.

<ul class="post-list">
  {% for post in site.posts %}
  <li>
    <span class="post-meta">{{ post.date | date: "%b %-d, %Y" }}</span>
    <a href="{{ post.url | relative_url }}">{{ post.title | escape }}</a>
  </li>
  {% endfor %}
</ul>

{% if site.posts.size == 0 %}
<p>No posts yet. Stay tuned!</p>
{% endif %}
