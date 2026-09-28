import sys
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np

# Set default font to Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

def pca_plot(df, group_names, group_colors, label_col=None):
    """
    Performs PCA on protein intensity data and creates a scatter plot.

    Args:
        df: DataFrame containing protein intensity data and group labels.
        group_names: List of group names (e.g., ["Normal", "HG", "HL"]).
        group_colors: List of colors corresponding to each group name.
        label_col: Column name with sample labels to display on the plot.

    Returns:
        None (displays the PCA plot).
    """

    # Separate features (protein intensities) and target (group)
    X = df.iloc[:, 2:]  # Select columns from the 3rd column onwards (protein intensities)
    y = df['Group']

    # Perform PCA
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

    # Calculate explained variance for axis labels
    explained_variance = pca.explained_variance_ratio_ * 100

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(5, 5))

    # Store text objects for interactive dragging
    text_objects = []
    
    for group, color in zip(group_names, group_colors):
        indicesToKeep = y == group
        ax.scatter(principalDf.loc[indicesToKeep, 'PC1'],
                   principalDf.loc[indicesToKeep, 'PC2'],
                   color=[color], label=group, s=50)

        # Add sample labels if specified
        if label_col:
            indices_list = list(principalDf[indicesToKeep].index)
            # Define 8 directions to spread labels around points
            directions = [
                (0, 1.2, 'center', 'bottom'),      # North
                (0, -1.2, 'center', 'top'),        # South
                (1.2, 0, 'left', 'center'),        # East
                (-1.2, 0, 'right', 'center'),      # West
                (0.85, 0.85, 'left', 'bottom'),    # Northeast
                (-0.85, 0.85, 'right', 'bottom'),  # Northwest
                (0.85, -0.85, 'left', 'top'),      # Southeast
                (-0.85, -0.85, 'right', 'top')     # Southwest
            ]
            
            for idx, i in enumerate(indices_list):
                # Cycle through directions to spread labels
                x_offset, y_offset, ha, va = directions[idx % len(directions)]
                
                text_obj = ax.text(principalDf.loc[i, 'PC1'] + x_offset, 
                principalDf.loc[i, 'PC2'] + y_offset,
                str(df.loc[i, label_col]), 
                fontsize=9, ha=ha, va=va, fontfamily='Arial', 
                picker=5)  # Enable picking for dragging (5 pixel tolerance)
                text_objects.append(text_obj)

    # Set axis labels with explained variance
    ax.set_xlabel(f'PC 1 ({explained_variance[0]:.2f}%)', fontsize=13, fontweight='bold', fontfamily='Arial')
    ax.set_ylabel(f'PC 2 ({explained_variance[1]:.2f}%)', fontsize=13, fontweight='bold', fontfamily='Arial')
    ax.set_title('PCA', fontsize=17, fontweight='bold', fontfamily='Arial')
    ax.legend(prop={'family': 'Arial'})
    
    # Interactive label dragging functionality
    if label_col and text_objects:
        selected_text = None
        click_offset = None
        
        def on_press(event):
            nonlocal selected_text, click_offset
            if event.inaxes != ax or event.button != 1:  # Only left mouse button
                return
            # Check if any text object was clicked
            for text_obj in text_objects:
                contains, _ = text_obj.contains(event)
                if contains:
                    selected_text = text_obj
                    # Get current text position
                    pos = text_obj.get_position()
                    # Calculate offset from click position to text position
                    click_offset = (event.xdata - pos[0], event.ydata - pos[1])
                    text_obj.set_alpha(0.7)  # Slightly fade when selected
                    fig.canvas.draw_idle()
                    break
        
        def on_motion(event):
            nonlocal selected_text, click_offset
            if selected_text is None or event.inaxes != ax:
                return
            if click_offset is not None:
                # Update text position based on mouse position minus offset
                new_x = event.xdata - click_offset[0]
                new_y = event.ydata - click_offset[1]
                selected_text.set_position((new_x, new_y))
                fig.canvas.draw_idle()
        
        def on_release(event):
            nonlocal selected_text, click_offset
            if selected_text is not None:
                selected_text.set_alpha(1.0)  # Restore full opacity
                fig.canvas.draw_idle()
            selected_text = None
            click_offset = None
        
        # Connect event handlers
        fig.canvas.mpl_connect('button_press_event', on_press)
        fig.canvas.mpl_connect('motion_notify_event', on_motion)
        fig.canvas.mpl_connect('button_release_event', on_release)
    
    plt.show()


if __name__ == "__main__":
    # Load data using your custom function
    df = pd.read_csv(r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\new_runs cell lysates 06-2025\04-PCA\PCA_new_run.csv")
    group_names = ["Normal", "HG", "HL","Combo"]
    group_colors = [(0,0, 1), (0.99, 0.57, 0), (0.99, 0.14, 0), (0.57, 0.125, 0.57)]
    # group_colors = ["cornflowerblue", "khaki","salmon","orchid"]
    sample_label_col = 'Sample'  # Replace with the column name for sample labels
    pca_plot(df, group_names, group_colors, label_col=sample_label_col)
