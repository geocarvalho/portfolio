---
layout: post
title: "Simple Text-to-Speech: Converting Text Files to Audio with Python"
date: 2026-01-21
categories: tools
tags: [python, text-to-speech, tts, cli, edge-tts, accessibility]
---

I often find myself wanting to listen to articles, notes, or documentation while commuting or exercising. Copy-pasting text into online TTS services felt clunky, so I built [simple-txt-to-speech](https://github.com/geocarvalho/simple-txt-to-speech) — a minimal command-line tool that converts any text file to MP3 using Microsoft Edge's neural TTS voices.

## Why Edge TTS?

Most TTS APIs are either paid (Google Cloud, AWS Polly) or produce robotic-sounding audio. Microsoft Edge TTS offers high-quality neural voices for free, with support for dozens of languages including Brazilian Portuguese. The [edge-tts](https://pypi.org/project/edge-tts/) Python package wraps this service cleanly, making it easy to integrate.

## How It Works

The tool reads a `.txt` file, sends the text to Edge TTS, and saves the resulting audio as `.mp3`. Under the hood it is straightforward:

```
Text File (.txt)
    │
    ├── Read & validate UTF-8 content
    │
    ├── Edge TTS API
    │     ├── Neural voice synthesis
    │     └── Streaming audio chunks
    │
    └── Output MP3 file
```

## Quick Start

Install and run in three steps:

```bash
git clone https://github.com/geocarvalho/simple-txt-to-speech.git
cd simple-txt-to-speech
pip install edge-tts
```

Basic conversion with default settings:

```bash
python bin/reader.py -i my_article.txt -o my_article.mp3
```

## Key Features

**Customizable speech rate** — speed up or slow down the output by a percentage. I usually set `+15%` for articles I already know well and `-10%` for dense technical content:

```bash
python bin/reader.py -i notes.txt -o notes.mp3 -r 15
```

**Multi-language voice selection** — list available voices filtered by locale. This is great for language learning or creating content in Portuguese:

```bash
python bin/reader.py --list-voices --locale pt-BR
```

Some voice options I use regularly:

| Voice | Language | Gender |
|---|---|---|
| `en-US-GuyNeural` | American English | Male |
| `en-US-AriaNeural` | American English | Female |
| `pt-BR-FranciscaNeural` | Brazilian Portuguese | Female |
| `en-GB-SoniaNeural` | British English | Female |

**Progress tracking** — a real-time progress bar shows conversion status, useful for longer texts:

```
[████████████████████░░░░░░░░░░░░░░░░░░] 60% (1.2 MB)
✓ Success! Audio saved to: result/result.mp3
```

## Use Cases

- **Commute listening** — convert blog posts or papers to audio for the drive
- **Accessibility** — make text content available as audio for visually impaired users
- **Language learning** — generate native-sounding audio from vocabulary lists or texts
- **Proofreading** — hearing your own writing read back helps catch awkward phrasing

## What's Next

Some ideas for future improvements:

- Support for PDF and Markdown input formats
- Chapter splitting for long texts
- Batch processing of multiple files
- A simple web interface

The tool is intentionally minimal — a single Python script with one dependency. Check out the [repository](https://github.com/geocarvalho/simple-txt-to-speech) for the full source and feel free to contribute.
