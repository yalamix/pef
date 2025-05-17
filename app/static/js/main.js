// Toggle mobile menu
document.addEventListener('DOMContentLoaded', function () {
  const toggleButton = document.getElementById('toggle-sidebar');
  const sidebar = document.getElementById('sidebar');

  if (toggleButton) {
    toggleButton.addEventListener('click', function () {
      sidebar.classList.toggle('hidden');
    });
  }

  // Handle responsive behavior
  function handleResize() {
    if (window.innerWidth < 768 && sidebar) {
      sidebar.classList.add('hidden');
    } else if (sidebar) {
      sidebar.classList.remove('hidden');
    }
  }

  // Initial check
  handleResize();

  // Listen for window resize
  window.addEventListener('resize', handleResize);
});

document.getElementById("openModal").addEventListener("click", function () {
  document.getElementById("modal").classList.remove("hidden");
});

document.getElementById("closeModal").addEventListener("click", function () {
  document.getElementById("modal").classList.add("hidden");
});

document.getElementById("cancelModal").addEventListener("click", function () {
  document.getElementById("modal").classList.add("hidden");
});

// Close modal when clicking outside of it
document.getElementById("modal").addEventListener("click", function (event) {
  if (event.target === this) {
    this.classList.add("hidden");
  }
});

// Dropdown portal logic
const dropdownButton = document.getElementById('dropdownButton');
const dropdownMenu = document.getElementById('dropdownMenu');
let dropdownAppended = false; // Track if menu is already appended

dropdownButton.addEventListener("click", function (e) {
  e.stopPropagation();
  if (dropdownMenu.classList.contains("hidden")) {
    // Calculate button's position
    const rect = dropdownButton.getBoundingClientRect();
    // Set dropdown menu styles for fixed positioning
    dropdownMenu.style.position = "fixed";
    dropdownMenu.style.top = rect.bottom + "px";
    dropdownMenu.style.left = rect.left + "px";
    // Append to body if not already
    if (!dropdownAppended) {
      document.body.appendChild(dropdownMenu);
      dropdownAppended = true;
    }
    dropdownMenu.classList.remove("hidden");
  } else {
    dropdownMenu.classList.add("hidden");
  }
});

// Update dropdown button text when an option is clicked, then close the dropdown
document.querySelectorAll('.dropdown-option').forEach(option => {
  option.addEventListener("click", function (e) {
    e.preventDefault();
    // Update button text with the clicked option's text
    dropdownSelectedText.textContent = this.textContent;
    // Close the dropdown
    dropdownMenu.classList.add("hidden");
  });
});

// Close dropdown if clicked outside
document.addEventListener("click", function (e) {
  if (!dropdownMenu.contains(e.target) && e.target !== dropdownButton) {
    dropdownMenu.classList.add("hidden");
  }
});