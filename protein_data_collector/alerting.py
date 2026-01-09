"""
Alerting system for the Protein Data Collector.

This module provides alerting capabilities for collection failures, API issues,
and system health problems with configurable notification channels.
"""

import asyncio
import logging
import smtplib
import json
from datetime import datetime, timedelta
from typing import Dict, Any, List, Optional, Callable
from dataclasses import dataclass, asdict
from enum import Enum
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import httpx

from .config import get_config
from .monitoring import HealthStatus, ComponentHealth, SystemHealth

logger = logging.getLogger(__name__)


class AlertSeverity(Enum):
    """Alert severity levels."""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


class AlertType(Enum):
    """Types of alerts."""
    COLLECTION_FAILURE = "collection_failure"
    API_FAILURE = "api_failure"
    DATABASE_FAILURE = "database_failure"
    SYSTEM_RESOURCE = "system_resource"
    DATA_QUALITY = "data_quality"
    HEALTH_CHECK = "health_check"


@dataclass
class Alert:
    """Alert data structure."""
    id: str
    type: AlertType
    severity: AlertSeverity
    title: str
    message: str
    details: Dict[str, Any]
    timestamp: datetime
    resolved: bool = False
    resolved_at: Optional[datetime] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert alert to dictionary."""
        return {
            "id": self.id,
            "type": self.type.value,
            "severity": self.severity.value,
            "title": self.title,
            "message": self.message,
            "details": self.details,
            "timestamp": self.timestamp.isoformat(),
            "resolved": self.resolved,
            "resolved_at": self.resolved_at.isoformat() if self.resolved_at else None
        }


class AlertChannel:
    """Base class for alert notification channels."""
    
    async def send_alert(self, alert: Alert) -> bool:
        """Send alert notification. Returns True if successful."""
        raise NotImplementedError


class LogAlertChannel(AlertChannel):
    """Alert channel that logs alerts."""
    
    def __init__(self):
        self.logger = logging.getLogger(f"{__name__}.alerts")
    
    async def send_alert(self, alert: Alert) -> bool:
        """Log the alert."""
        try:
            log_level = {
                AlertSeverity.LOW: logging.INFO,
                AlertSeverity.MEDIUM: logging.WARNING,
                AlertSeverity.HIGH: logging.ERROR,
                AlertSeverity.CRITICAL: logging.CRITICAL
            }.get(alert.severity, logging.WARNING)
            
            self.logger.log(
                log_level,
                f"ALERT: {alert.title}",
                extra={
                    "alert_id": alert.id,
                    "alert_type": alert.type.value,
                    "alert_severity": alert.severity.value,
                    "alert_message": alert.message,
                    "alert_details": alert.details,
                    "alert_timestamp": alert.timestamp.isoformat()
                }
            )
            return True
        except Exception as e:
            logger.error(f"Failed to log alert: {e}")
            return False


class EmailAlertChannel(AlertChannel):
    """Alert channel that sends email notifications."""
    
    def __init__(self, smtp_host: str, smtp_port: int, username: str, 
                 password: str, from_email: str, to_emails: List[str],
                 use_tls: bool = True):
        self.smtp_host = smtp_host
        self.smtp_port = smtp_port
        self.username = username
        self.password = password
        self.from_email = from_email
        self.to_emails = to_emails
        self.use_tls = use_tls
    
    async def send_alert(self, alert: Alert) -> bool:
        """Send alert via email."""
        try:
            # Create message
            msg = MIMEMultipart()
            msg['From'] = self.from_email
            msg['To'] = ', '.join(self.to_emails)
            msg['Subject'] = f"[{alert.severity.value.upper()}] {alert.title}"
            
            # Create email body
            body = self._create_email_body(alert)
            msg.attach(MIMEText(body, 'html'))
            
            # Send email
            with smtplib.SMTP(self.smtp_host, self.smtp_port) as server:
                if self.use_tls:
                    server.starttls()
                server.login(self.username, self.password)
                server.send_message(msg)
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to send email alert: {e}")
            return False
    
    def _create_email_body(self, alert: Alert) -> str:
        """Create HTML email body."""
        severity_colors = {
            AlertSeverity.LOW: "#28a745",
            AlertSeverity.MEDIUM: "#ffc107",
            AlertSeverity.HIGH: "#fd7e14",
            AlertSeverity.CRITICAL: "#dc3545"
        }
        
        color = severity_colors.get(alert.severity, "#6c757d")
        
        return f"""
        <html>
        <body>
            <h2 style="color: {color};">Protein Data Collector Alert</h2>
            <p><strong>Severity:</strong> <span style="color: {color};">{alert.severity.value.upper()}</span></p>
            <p><strong>Type:</strong> {alert.type.value.replace('_', ' ').title()}</p>
            <p><strong>Time:</strong> {alert.timestamp.strftime('%Y-%m-%d %H:%M:%S UTC')}</p>
            <p><strong>Message:</strong> {alert.message}</p>
            
            <h3>Details:</h3>
            <pre style="background-color: #f8f9fa; padding: 10px; border-radius: 5px;">
{json.dumps(alert.details, indent=2)}
            </pre>
            
            <hr>
            <p><small>This alert was generated by the Protein Data Collector monitoring system.</small></p>
        </body>
        </html>
        """


class WebhookAlertChannel(AlertChannel):
    """Alert channel that sends webhooks (e.g., to Slack, Discord, etc.)."""
    
    def __init__(self, webhook_url: str, headers: Optional[Dict[str, str]] = None):
        self.webhook_url = webhook_url
        self.headers = headers or {}
    
    async def send_alert(self, alert: Alert) -> bool:
        """Send alert via webhook."""
        try:
            payload = {
                "text": f"🚨 {alert.title}",
                "attachments": [
                    {
                        "color": self._get_color(alert.severity),
                        "fields": [
                            {"title": "Severity", "value": alert.severity.value.upper(), "short": True},
                            {"title": "Type", "value": alert.type.value.replace('_', ' ').title(), "short": True},
                            {"title": "Time", "value": alert.timestamp.strftime('%Y-%m-%d %H:%M:%S UTC'), "short": True},
                            {"title": "Message", "value": alert.message, "short": False}
                        ]
                    }
                ]
            }
            
            async with httpx.AsyncClient() as client:
                response = await client.post(
                    self.webhook_url,
                    json=payload,
                    headers=self.headers,
                    timeout=10.0
                )
                response.raise_for_status()
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to send webhook alert: {e}")
            return False
    
    def _get_color(self, severity: AlertSeverity) -> str:
        """Get color for severity level."""
        colors = {
            AlertSeverity.LOW: "good",
            AlertSeverity.MEDIUM: "warning",
            AlertSeverity.HIGH: "danger",
            AlertSeverity.CRITICAL: "danger"
        }
        return colors.get(severity, "warning")


class AlertManager:
    """Manages alerts and notifications."""
    
    def __init__(self):
        self.channels: List[AlertChannel] = []
        self.active_alerts: Dict[str, Alert] = {}
        self.alert_history: List[Alert] = []
        self.alert_rules: List[Callable[[Dict[str, Any]], Optional[Alert]]] = []
        self.config = get_config()
        
        # Initialize default channels
        self._initialize_channels()
        self._initialize_rules()
    
    def _initialize_channels(self):
        """Initialize alert channels based on configuration."""
        # Always add log channel
        self.channels.append(LogAlertChannel())
        
        # Add email channel if configured
        # This would be configured via environment variables or config file
        # For now, we'll just add the log channel
    
    def _initialize_rules(self):
        """Initialize alert rules."""
        self.alert_rules = [
            self._check_collection_failures,
            self._check_api_failures,
            self._check_database_health,
            self._check_system_resources,
            self._check_data_quality
        ]
    
    async def process_system_health(self, health: SystemHealth):
        """Process system health and generate alerts if needed."""
        for component in health.components:
            if component.status == HealthStatus.UNHEALTHY:
                await self._create_alert(
                    AlertType.HEALTH_CHECK,
                    AlertSeverity.HIGH,
                    f"{component.name.title()} Component Unhealthy",
                    f"Component {component.name} is unhealthy: {component.error_message}",
                    {
                        "component": component.name,
                        "status": component.status.value,
                        "error": component.error_message,
                        "response_time_ms": component.response_time_ms,
                        "details": component.details
                    }
                )
            elif component.status == HealthStatus.DEGRADED:
                await self._create_alert(
                    AlertType.HEALTH_CHECK,
                    AlertSeverity.MEDIUM,
                    f"{component.name.title()} Component Degraded",
                    f"Component {component.name} is experiencing degraded performance",
                    {
                        "component": component.name,
                        "status": component.status.value,
                        "response_time_ms": component.response_time_ms,
                        "details": component.details
                    }
                )
    
    async def process_collection_failure(self, phase: str, error: str, details: Dict[str, Any]):
        """Process collection failure and generate alert."""
        await self._create_alert(
            AlertType.COLLECTION_FAILURE,
            AlertSeverity.HIGH,
            f"Data Collection Failed: {phase}",
            f"Data collection failed during {phase}: {error}",
            {
                "phase": phase,
                "error": error,
                **details
            }
        )
    
    async def process_api_failure(self, api_name: str, endpoint: str, error: str, 
                                 failure_count: int, details: Dict[str, Any]):
        """Process API failure and generate alert."""
        severity = AlertSeverity.CRITICAL if failure_count > 10 else AlertSeverity.HIGH
        
        await self._create_alert(
            AlertType.API_FAILURE,
            severity,
            f"API Failure: {api_name}",
            f"Multiple failures detected for {api_name} API at {endpoint}. "
            f"Failure count: {failure_count}. Error: {error}",
            {
                "api_name": api_name,
                "endpoint": endpoint,
                "error": error,
                "failure_count": failure_count,
                **details
            }
        )
    
    async def process_data_quality_issue(self, issue_type: str, count: int, 
                                       total: int, details: Dict[str, Any]):
        """Process data quality issue and generate alert."""
        percentage = (count / total * 100) if total > 0 else 0
        severity = AlertSeverity.HIGH if percentage > 10 else AlertSeverity.MEDIUM
        
        await self._create_alert(
            AlertType.DATA_QUALITY,
            severity,
            f"Data Quality Issue: {issue_type}",
            f"Data quality issue detected: {count}/{total} records ({percentage:.1f}%) "
            f"affected by {issue_type}",
            {
                "issue_type": issue_type,
                "affected_count": count,
                "total_count": total,
                "percentage": percentage,
                **details
            }
        )
    
    async def _create_alert(self, alert_type: AlertType, severity: AlertSeverity,
                           title: str, message: str, details: Dict[str, Any]):
        """Create and send an alert."""
        alert_id = f"{alert_type.value}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        alert = Alert(
            id=alert_id,
            type=alert_type,
            severity=severity,
            title=title,
            message=message,
            details=details,
            timestamp=datetime.now()
        )
        
        # Store alert
        self.active_alerts[alert_id] = alert
        self.alert_history.append(alert)
        
        # Send notifications
        await self._send_alert(alert)
    
    async def _send_alert(self, alert: Alert):
        """Send alert through all configured channels."""
        tasks = []
        for channel in self.channels:
            tasks.append(channel.send_alert(alert))
        
        if tasks:
            results = await asyncio.gather(*tasks, return_exceptions=True)
            success_count = sum(1 for result in results if result is True)
            logger.info(f"Alert {alert.id} sent to {success_count}/{len(tasks)} channels")
    
    async def resolve_alert(self, alert_id: str):
        """Mark an alert as resolved."""
        if alert_id in self.active_alerts:
            alert = self.active_alerts[alert_id]
            alert.resolved = True
            alert.resolved_at = datetime.now()
            del self.active_alerts[alert_id]
            logger.info(f"Alert {alert_id} resolved")
    
    def get_active_alerts(self) -> List[Alert]:
        """Get all active alerts."""
        return list(self.active_alerts.values())
    
    def get_alert_history(self, hours: int = 24) -> List[Alert]:
        """Get alert history for the specified number of hours."""
        cutoff = datetime.now() - timedelta(hours=hours)
        return [alert for alert in self.alert_history if alert.timestamp >= cutoff]
    
    def _check_collection_failures(self, metrics: Dict[str, Any]) -> Optional[Alert]:
        """Check for collection failures."""
        # This would be implemented based on actual metrics
        return None
    
    def _check_api_failures(self, metrics: Dict[str, Any]) -> Optional[Alert]:
        """Check for API failures."""
        # This would be implemented based on actual metrics
        return None
    
    def _check_database_health(self, metrics: Dict[str, Any]) -> Optional[Alert]:
        """Check database health."""
        # This would be implemented based on actual metrics
        return None
    
    def _check_system_resources(self, metrics: Dict[str, Any]) -> Optional[Alert]:
        """Check system resource usage."""
        # This would be implemented based on actual metrics
        return None
    
    def _check_data_quality(self, metrics: Dict[str, Any]) -> Optional[Alert]:
        """Check data quality metrics."""
        # This would be implemented based on actual metrics
        return None


# Global alert manager instance
_alert_manager = None


def get_alert_manager() -> AlertManager:
    """Get the global alert manager instance."""
    global _alert_manager
    if _alert_manager is None:
        _alert_manager = AlertManager()
    return _alert_manager