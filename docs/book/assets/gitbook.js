(function () {
  "use strict";

  var bookEl = document.querySelector(".book");
  var toggleBtn = document.getElementById("sidebar-toggle");
  var searchInput = document.getElementById("book-search-input");
  var searchToggleBtn = document.getElementById("search-toggle-btn");
  var searchResultsEl = document.getElementById("search-results");
  var fontBtn = document.getElementById("font-settings-btn");
  var fontMenu = document.getElementById("font-settings-menu");
  var searchIndex = null;

  // Restore saved preferences
  var savedTheme = localStorage.getItem("chpy_theme") || "white";
  var savedFont = localStorage.getItem("chpy_font") || "sans";
  var savedSize = parseInt(localStorage.getItem("chpy_size") || "16", 10);
  var savedSidebar = localStorage.getItem("chpy_sidebar");

  function applyTheme(theme) {
    document.body.classList.remove("theme-white", "theme-sepia", "theme-night");
    document.body.classList.add("theme-" + theme);
    localStorage.setItem("chpy_theme", theme);
    document.querySelectorAll("[data-theme]").forEach(function (b) {
      b.classList.toggle("active", b.getAttribute("data-theme") === theme);
    });
  }

  function applyFont(font) {
    document.body.classList.toggle("font-serif", font === "serif");
    localStorage.setItem("chpy_font", font);
    document.querySelectorAll("[data-font]").forEach(function (b) {
      b.classList.toggle("active", b.getAttribute("data-font") === font);
    });
  }

  function applySize(sz) {
    var clamped = Math.max(13, Math.min(22, sz));
    document.documentElement.style.setProperty("--font-size-base", clamped + "px");
    localStorage.setItem("chpy_size", String(clamped));
  }

  applyTheme(savedTheme);
  applyFont(savedFont);
  applySize(savedSize);

  if (savedSidebar === "closed" && window.innerWidth > 960 && bookEl) {
    bookEl.classList.add("without-summary");
  }

  if (toggleBtn && bookEl) {
    toggleBtn.addEventListener("click", function () {
      if (window.innerWidth <= 960) {
        bookEl.classList.toggle("mobile-summary-open");
      } else {
        bookEl.classList.toggle("without-summary");
        localStorage.setItem("chpy_sidebar", bookEl.classList.contains("without-summary") ? "closed" : "open");
      }
    });
  }

  // Font settings dropdown
  if (fontBtn && fontMenu) {
    fontBtn.addEventListener("click", function (e) {
      e.stopPropagation();
      fontMenu.classList.toggle("open");
    });
    document.addEventListener("click", function (e) {
      if (!fontMenu.contains(e.target) && e.target !== fontBtn) {
        fontMenu.classList.remove("open");
      }
    });
  }

  document.querySelectorAll("[data-theme]").forEach(function (btn) {
    btn.addEventListener("click", function () { applyTheme(btn.getAttribute("data-theme")); });
  });
  document.querySelectorAll("[data-font]").forEach(function (btn) {
    btn.addEventListener("click", function () { applyFont(btn.getAttribute("data-font")); });
  });
  var decBtn = document.getElementById("font-size-dec");
  var incBtn = document.getElementById("font-size-inc");
  if (decBtn) decBtn.addEventListener("click", function () {
    applySize(parseInt(localStorage.getItem("chpy_size") || "16", 10) - 1);
  });
  if (incBtn) incBtn.addEventListener("click", function () {
    applySize(parseInt(localStorage.getItem("chpy_size") || "16", 10) + 1);
  });

  // Search
  function loadSearchIndex(cb) {
    if (searchIndex) return cb(searchIndex);
    fetch("assets/search_index.json")
      .then(function (r) { return r.json(); })
      .then(function (data) { searchIndex = data; cb(data); })
      .catch(function () { cb([]); });
  }

  function escapeHtml(s) {
    return s.replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;", "'": "&#39;" }[c];
    });
  }

  function runSearch(query) {
    var q = (query || "").trim().toLowerCase();
    if (!searchResultsEl) return;
    if (q.length < 2) {
      searchResultsEl.classList.remove("active");
      searchResultsEl.innerHTML = "";
      return;
    }
    loadSearchIndex(function (items) {
      var terms = q.split(/\s+/).filter(Boolean);
      var matches = [];
      items.forEach(function (item) {
        var hay = (item.title + " " + item.text).toLowerCase();
        if (terms.every(function (t) { return hay.indexOf(t) !== -1; })) {
          var idx = item.text.toLowerCase().indexOf(terms[0]);
          var start = Math.max(0, idx - 60);
          var snippet = (start > 0 ? "..." : "") + item.text.slice(start, start + 180) + "...";
          matches.push({ item: item, snippet: snippet });
        }
      });
      var html = "<strong>Search results (" + matches.length + ")</strong>";
      if (matches.length === 0) {
        html += "<p style='margin:6px 0 0;color:var(--text-muted)'>No matching sections found.</p>";
      } else {
        matches.slice(0, 20).forEach(function (m) {
          html += "<div class='search-result-item'><a href='" + escapeHtml(m.item.url) + "'>" +
            escapeHtml(m.item.title) + "</a><div class='search-result-snippet'>" +
            escapeHtml(m.snippet) + "</div></div>";
        });
      }
      searchResultsEl.innerHTML = html;
      searchResultsEl.classList.add("active");
    });
  }

  if (searchInput) {
    searchInput.addEventListener("input", function () { runSearch(searchInput.value); });
  }
  if (searchToggleBtn && searchInput) {
    searchToggleBtn.addEventListener("click", function () {
      if (bookEl && window.innerWidth <= 960) {
        bookEl.classList.add("mobile-summary-open");
      } else if (bookEl && bookEl.classList.contains("without-summary")) {
        bookEl.classList.remove("without-summary");
      }
      searchInput.focus();
    });
  }

  // Copy buttons for code blocks
  document.querySelectorAll(".code-wrapper pre").forEach(function (pre) {
    var wrapper = pre.parentElement;
    var btn = document.createElement("button");
    btn.className = "copy-code-btn";
    btn.type = "button";
    btn.textContent = "Copy";
    btn.addEventListener("click", function () {
      var codeText = pre.innerText || pre.textContent || "";
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(codeText).then(function () {
          btn.textContent = "Copied!";
          setTimeout(function () { btn.textContent = "Copy"; }, 1500);
        });
      }
    });
    wrapper.appendChild(btn);
  });

  // Keyboard navigation (Left/Right arrows outside inputs & modals)
  document.addEventListener("keydown", function (e) {
    if (e.target && (e.target.tagName === "INPUT" || e.target.tagName === "TEXTAREA")) return;
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    var lb = document.getElementById("gallery-lightbox");
    if (lb && !lb.hidden) return;
    if (e.key === "ArrowLeft") {
      var prev = document.querySelector(".navigation-prev");
      if (prev && prev.getAttribute("href")) window.location.href = prev.getAttribute("href");
    } else if (e.key === "ArrowRight") {
      var next = document.querySelector(".navigation-next");
      if (next && next.getAttribute("href")) window.location.href = next.getAttribute("href");
    }
  });

  // Figure Gallery filtering, thumbnail density, and lightbox
  function initGallery() {
    var root = document.getElementById("gallery-root");
    if (!root) return;

    var filterInput = document.getElementById("gallery-filter-input");
    var countBadge = document.getElementById("gallery-count-badge");
    var emptyMsg = document.getElementById("gallery-empty-msg");
    var pills = Array.prototype.slice.call(document.querySelectorAll(".gallery-pill"));
    var sizeBtns = Array.prototype.slice.call(document.querySelectorAll(".gallery-size-btn"));
    var groups = Array.prototype.slice.call(document.querySelectorAll(".gallery-group"));
    var cards = Array.prototype.slice.call(document.querySelectorAll(".gallery-card"));

    var lightbox = document.getElementById("gallery-lightbox");
    var lbBackdrop = lightbox ? lightbox.querySelector(".lightbox-backdrop") : null;
    var lbClose = document.getElementById("lightbox-close");
    var lbPrev = document.getElementById("lightbox-prev");
    var lbNext = document.getElementById("lightbox-next");
    var lbImg = document.getElementById("lightbox-img");
    var lbCounter = document.getElementById("lightbox-counter");
    var lbCaption = document.getElementById("lightbox-caption");
    var lbMeta = document.getElementById("lightbox-meta");
    var lbJump = document.getElementById("lightbox-jump");
    var lbRaw = document.getElementById("lightbox-raw");

    var activeCh = "all";
    var totalCount = cards.length;
    var visibleCards = cards.slice();
    var currentLbPos = -1;

    function applyGridSize(sz) {
      var mode = (sz === "compact" || sz === "large") ? sz : "medium";
      root.classList.remove("gallery-grid-compact", "gallery-grid-large");
      if (mode !== "medium") {
        root.classList.add("gallery-grid-" + mode);
      }
      sizeBtns.forEach(function (b) {
        b.classList.toggle("active", b.getAttribute("data-grid-size") === mode);
      });
      localStorage.setItem("chpy_gallery_size", mode);
    }

    applyGridSize(localStorage.getItem("chpy_gallery_size") || "medium");
    sizeBtns.forEach(function (btn) {
      btn.addEventListener("click", function () {
        applyGridSize(btn.getAttribute("data-grid-size"));
      });
    });

    function applyFilters() {
      var rawQ = filterInput ? filterInput.value.trim().toLowerCase() : "";
      var terms = rawQ ? rawQ.split(/\s+/).filter(Boolean) : [];
      var shown = 0;
      visibleCards = [];

      groups.forEach(function (grp) {
        var grpCh = grp.getAttribute("data-ch-group");
        var grpCards = Array.prototype.slice.call(grp.querySelectorAll(".gallery-card"));
        var grpShown = 0;

        grpCards.forEach(function (card) {
          var chMatch = activeCh === "all" || grpCh === activeCh;
          var textMatch = true;
          if (terms.length > 0) {
            var hay = card.getAttribute("data-search") || "";
            textMatch = terms.every(function (t) { return hay.indexOf(t) !== -1; });
          }
          var vis = chMatch && textMatch;
          card.hidden = !vis;
          if (vis) {
            grpShown += 1;
            shown += 1;
            visibleCards.push(card);
          }
        });
        grp.hidden = grpShown === 0;
      });

      if (countBadge) {
        countBadge.textContent = shown + " / " + totalCount + " figures";
      }
      if (emptyMsg) {
        emptyMsg.hidden = shown > 0;
      }
    }

    function setActiveChapter(chKey) {
      activeCh = chKey || "all";
      pills.forEach(function (p) {
        p.classList.toggle("active", p.getAttribute("data-ch-filter") === activeCh);
      });
      applyFilters();
    }

    pills.forEach(function (pill) {
      pill.addEventListener("click", function () {
        setActiveChapter(pill.getAttribute("data-ch-filter"));
      });
    });

    if (filterInput) {
      filterInput.addEventListener("input", applyFilters);
    }

    window.addEventListener("hashchange", function () {
      if (window.location.hash && window.location.hash.indexOf("#gallery-ch-") === 0) {
        var target = document.getElementById(window.location.hash.slice(1));
        if (target && target.closest(".gallery-group") && target.closest(".gallery-group").hidden) {
          setActiveChapter("all");
          target.scrollIntoView();
        }
      }
    });

    function showLightboxAt(pos) {
      if (!lightbox || visibleCards.length === 0) return;
      if (pos < 0) pos = visibleCards.length - 1;
      if (pos >= visibleCards.length) pos = 0;
      currentLbPos = pos;

      var btn = visibleCards[pos].querySelector(".gallery-thumb-btn");
      if (!btn) return;
      var src = btn.getAttribute("data-full-src") || "";
      var caption = btn.getAttribute("data-caption") || "";
      var chapter = btn.getAttribute("data-chapter") || "";
      var section = btn.getAttribute("data-section") || "";
      var url = btn.getAttribute("data-url") || "#";
      var fname = btn.getAttribute("data-file") || "";

      if (lbImg) {
        lbImg.src = src;
        lbImg.alt = caption;
      }
      if (lbCounter) lbCounter.textContent = "Figure " + (pos + 1) + " of " + visibleCards.length;
      if (lbCaption) lbCaption.textContent = caption;
      if (lbMeta) lbMeta.textContent = chapter + " \u2022 " + section + " (" + fname + ")";
      if (lbJump) lbJump.href = url;
      if (lbRaw) lbRaw.href = src;

      lightbox.hidden = false;
      document.body.style.overflow = "hidden";
    }

    function closeLightbox() {
      if (!lightbox) return;
      lightbox.hidden = true;
      document.body.style.overflow = "";
    }

    cards.forEach(function (card) {
      var btn = card.querySelector(".gallery-thumb-btn");
      if (!btn) return;
      btn.addEventListener("click", function () {
        var idx = visibleCards.indexOf(card);
        if (idx === -1) {
          applyFilters();
          idx = visibleCards.indexOf(card);
        }
        if (idx !== -1) showLightboxAt(idx);
      });
    });

    if (lbBackdrop) lbBackdrop.addEventListener("click", closeLightbox);
    if (lbClose) lbClose.addEventListener("click", closeLightbox);
    if (lbPrev) lbPrev.addEventListener("click", function () { showLightboxAt(currentLbPos - 1); });
    if (lbNext) lbNext.addEventListener("click", function () { showLightboxAt(currentLbPos + 1); });

    document.addEventListener("keydown", function (e) {
      if (!lightbox || lightbox.hidden) return;
      if (e.key === "Escape") {
        e.preventDefault();
        closeLightbox();
      } else if (e.key === "ArrowLeft") {
        e.preventDefault();
        showLightboxAt(currentLbPos - 1);
      } else if (e.key === "ArrowRight") {
        e.preventDefault();
        showLightboxAt(currentLbPos + 1);
      }
    });
  }

  initGallery();
})();
